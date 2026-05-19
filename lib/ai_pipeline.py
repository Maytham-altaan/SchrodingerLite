"""AI Discovery pipeline — closed-loop de novo drug design.

    Round 1:  Claude designs N novel molecules for the target
              → filter (validity, drug-likeness, SA, PAINS, novelty)
              → ligprep (3D + PDBQT)
              → glide dock against the prepared receptor
              → rank by ΔG

    Round 2..K:  Claude is shown the top hits + their scores
                 → designs N more (analogues of winners, plus new scaffolds)
                 → same filter / dock / rank pipeline
                 → merged ranked list

    Output:   ranked CSV of every molecule across all rounds,
              SDF poses of the top K winners,
              Claude's med-chem summary of the final top hits.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

from .utils import banner, get_logger, new_project
from . import ai_design, ai_filter, ligprep, docking, protprep


# ---------------------------------------------------------------------------

def _save_smi(molecules: list[dict], path: Path):
    with open(path, "w") as f:
        for m in molecules:
            name = m.get("name") or m["smiles"][:10]
            f.write(f"{m['smiles']}\t{name}\n")


def _save_jsonl(molecules: list[dict], path: Path):
    with open(path, "w") as f:
        for m in molecules:
            f.write(json.dumps(m) + "\n")


def _load_dock_csv(csv_path: Path) -> dict[str, float]:
    """Map ligand_name → best_score from a docking_results.csv."""
    out = {}
    if not csv_path.exists():
        return out
    for row in csv.DictReader(csv_path.open()):
        try:
            out[row["ligand"]] = float(row["best_score"])
        except (KeyError, ValueError, TypeError):
            pass
    return out


# ---------------------------------------------------------------------------

def run_discovery(target_description: str,
                  receptor_pdbqt: Path,
                  output_dir: Path,
                  ref_ligand: Path | None = None,
                  reference_library: Path | None = None,
                  rounds: int = 3,
                  per_round: int = 20,
                  top_k_feedback: int = 5,
                  exhaustiveness: int = 8,
                  grid_center: tuple[float,float,float] | None = None,
                  grid_size: tuple[float,float,float] = (22, 22, 22),
                  grid_padding: float = 8.0,
                  sa_max: float = 6.0,
                  novelty_min: float = 0.4) -> Path:
    """Run the full AI discovery loop. Returns final ranked CSV path."""
    log = get_logger("aidiscover")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "rounds").mkdir(exist_ok=True)

    # Grid box
    if ref_ligand and ref_ligand.exists():
        center, size = docking.grid_from_ligand(ref_ligand, padding=grid_padding)
        log.info(f"Grid box from ref ligand: center={center}, size={size}")
    elif grid_center is not None:
        center, size = grid_center, grid_size
        log.info(f"Grid box explicit: center={center}, size={size}")
    else:
        raise ValueError("Provide either ref_ligand or grid_center.")

    # Reference fingerprints for novelty
    ref_fps = ai_filter.load_reference_fps(reference_library) if reference_library else []
    if ref_fps:
        log.info(f"Loaded {len(ref_fps)} reference fingerprints for novelty scoring")

    all_molecules: list[dict] = []   # accumulated across rounds
    seen_canonical: set[str] = set() # canonical SMILES already kept
    used_names: set[str] = set()     # globally-unique ligand names
    feedback: list[dict] | None = None

    for r in range(1, rounds + 1):
        rdir = output_dir / "rounds" / f"round_{r:02d}"
        rdir.mkdir(parents=True, exist_ok=True)
        log.info(f"━━━━━━━━━━ ROUND {r}/{rounds} ━━━━━━━━━━")

        # 1. Generate with Claude
        proposed = ai_design.generate_round(
            target_description=target_description,
            n_molecules=per_round,
            feedback=feedback,
            round_num=r,
        )
        _save_jsonl(proposed, rdir / "01_proposed.jsonl")

        if not proposed:
            log.warning("Claude returned no parsable molecules this round; skipping.")
            continue

        # 2. Filter (canonicalizes SMILES into m["smiles"])
        log.info(f"Filtering {len(proposed)} proposals…")
        kept = ai_filter.filter_batch(proposed, ref_fps=ref_fps,
                                      sa_max=sa_max,
                                      novelty_min=novelty_min)
        log.info(f"  {len(kept)}/{len(proposed)} passed all filters")

        # 2b. Cross-round dedup by canonical SMILES + globally-unique names
        deduped: list[dict] = []
        dup_smiles = 0
        renamed = 0
        for m in kept:
            csmi = m["smiles"]
            if csmi in seen_canonical:
                dup_smiles += 1
                continue
            seen_canonical.add(csmi)
            name = m.get("name") or csmi[:10]
            if name in used_names:
                base = name
                k = 2
                while f"{base}_r{r}-{k}" in used_names:
                    k += 1
                name = f"{base}_r{r}-{k}"
                renamed += 1
            used_names.add(name)
            m["name"] = name
            deduped.append(m)
        if dup_smiles or renamed:
            log.info(f"  dedup: dropped {dup_smiles} duplicate SMILES, "
                     f"renamed {renamed} name collisions")
        kept = deduped
        _save_jsonl(kept, rdir / "02_filtered.jsonl")

        if not kept:
            log.warning("No molecules survived filters this round.")
            continue

        # 3. LigPrep
        smi_file = rdir / "ligands.smi"
        _save_smi(kept, smi_file)
        prep_dir = rdir / "ligprep"
        log.info("Preparing 3D ligands + PDBQT (LigPrep)…")
        ligprep.prepare_ligands(smi_file, prep_dir,
                                n_confs=8, make_pdbqt=True)
        pdbqts = list(prep_dir.glob("*.pdbqt"))
        if not pdbqts:
            log.warning("LigPrep produced no PDBQTs; skipping round.")
            continue

        # 4. Glide / Vina dock
        dock_dir = rdir / "docking"
        log.info(f"Docking {len(pdbqts)} ligands with AutoDock Vina…")
        csv_path = docking.dock_batch(
            receptor_pdbqt=receptor_pdbqt,
            ligand_files=pdbqts,
            output_dir=dock_dir,
            center=center, size=size,
            backend="vina",
            exhaustiveness=exhaustiveness,
            n_poses=5,
        )

        # 5. Merge scores back onto molecule records
        scores = _load_dock_csv(csv_path)
        for m in kept:
            name = m.get("name") or m["smiles"][:10]
            m["score"] = scores.get(name)
            m["round"] = r
            all_molecules.append(m)

        # 6. Build feedback for next round: top K hits
        scored = [m for m in kept if m.get("score") is not None]
        scored.sort(key=lambda m: m["score"])
        top = scored[:top_k_feedback]
        feedback = [
            {"smiles": m["smiles"],
             "score": round(m["score"], 2),
             "notes": (m.get("rationale","")[:80])}
            for m in top
        ]
        for t in top:
            log.info(f"  ★ {t.get('name','?'):20s}  {t['score']:.2f} kcal/mol  "
                     f"SA={t.get('SA','?')}  novelty={t.get('novelty','?')}")

    # ---- Final ranking across all rounds ---------------------------------

    scored_all = [m for m in all_molecules if m.get("score") is not None]
    scored_all.sort(key=lambda m: m["score"])

    final_csv = output_dir / "final_ranked.csv"
    fieldnames = ["rank", "name", "smiles", "score", "round",
                  "MW", "LogP", "HBD", "HBA", "TPSA", "SA",
                  "QED", "novelty", "rationale"]
    with open(final_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for i, m in enumerate(scored_all, 1):
            w.writerow({"rank": i, **m})
    log.info(f"Final ranked CSV → {final_csv}")

    # ---- Claude med-chem summary -----------------------------------------
    if scored_all:
        try:
            summary = ai_design.explain_results(scored_all[:10], target_description)
            (output_dir / "summary.txt").write_text(summary)
            log.info("Med-chem summary written → summary.txt")
        except Exception as e:
            log.warning(f"Summary generation failed: {e}")

    return final_csv


# ---------------------------------------------------------------------------
# CLI

def main(argv=None):
    p = argparse.ArgumentParser(
        prog="aidiscover",
        description="AI-driven de novo drug discovery (Claude + RDKit + Vina)."
    )
    p.add_argument("--target", required=True,
                   help="Plain-text description of the target & binding site, "
                        "OR path to a .txt file containing the description.")
    p.add_argument("-r", "--receptor", required=True, type=Path,
                   help="Prepared receptor PDBQT (from prepwizard)")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--ref-ligand", type=Path, help="Reference ligand PDB for grid box")
    g.add_argument("--center", nargs=3, type=float, metavar=("X","Y","Z"))
    p.add_argument("--size", nargs=3, type=float, default=[22,22,22],
                   metavar=("X","Y","Z"))
    p.add_argument("--padding", type=float, default=8.0)
    p.add_argument("-o", "--output", type=Path, default=None)
    p.add_argument("--rounds", type=int, default=3,
                   help="Discovery iterations (default 3)")
    p.add_argument("--per-round", type=int, default=20,
                   help="Molecules generated per round (default 20)")
    p.add_argument("--top-k", type=int, default=5,
                   help="Top hits fed back to Claude (default 5)")
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--reference-library", type=Path, default=None,
                   help="SMILES file of known actives — used for novelty scoring")
    p.add_argument("--sa-max", type=float, default=6.0)
    p.add_argument("--novelty-min", type=float, default=0.4)

    args = p.parse_args(argv)

    # Resolve target text
    tgt = Path(args.target)
    target_description = tgt.read_text() if tgt.exists() else args.target

    banner("AI Drug Discovery — Claude × Vina",
           f"rounds={args.rounds}  per_round={args.per_round}  "
           f"top_k={args.top_k}")
    print(f"  Target:\n    {target_description[:200]}{'…' if len(target_description)>200 else ''}")

    out = (new_project("aidiscover") / "output") if args.output is None else args.output
    out.mkdir(parents=True, exist_ok=True)

    final = run_discovery(
        target_description=target_description,
        receptor_pdbqt=args.receptor,
        output_dir=out,
        ref_ligand=args.ref_ligand,
        reference_library=args.reference_library,
        rounds=args.rounds,
        per_round=args.per_round,
        top_k_feedback=args.top_k,
        exhaustiveness=args.exhaustiveness,
        grid_center=tuple(args.center) if args.center else None,
        grid_size=tuple(args.size),
        grid_padding=args.padding,
        sa_max=args.sa_max,
        novelty_min=args.novelty_min,
    )

    print(f"\n✓ Discovery complete.")
    print(f"  Ranked candidates: {final}")
    print(f"  Med-chem summary:  {final.parent / 'summary.txt'}")
    print(f"  Per-round results: {final.parent / 'rounds/'}")


if __name__ == "__main__":
    main()
