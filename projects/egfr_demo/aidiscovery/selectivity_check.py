"""Wild-type selectivity check for the AI-discovered EGFR T790M/L858R hits.

For the top-N unique hits from the mutant-EGFR campaign:
  1. Prepare wild-type EGFR (PDB 1M17 — erlotinib co-crystal).
  2. Re-dock each ligand into the wild-type ATP-binding site at the same
     exhaustiveness as the mutant campaign.
  3. Compute ΔΔG = wt_score − mut_score.   Positive ΔΔG = mutant-selective.
  4. Write selectivity.csv and a side-by-side 2-pane HTML viewer.

Run:
    python selectivity_check.py --top 5
"""
from __future__ import annotations
import argparse
import base64
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]                       # .../SchrodingerLite
sys.path.insert(0, str(ROOT))
from lib import protprep, docking            # noqa: E402
from lib.viewer import _ensure_3dmol_js      # noqa: E402

# --- Paths -----------------------------------------------------------------
MUT_RECEPTOR    = ROOT / "projects" / "egfr_demo" / "receptor" / "protein_prepared.pdb"
MUT_RECEPTOR_QT = ROOT / "projects" / "egfr_demo" / "receptor" / "receptor.pdbqt"
ROUNDS_DIR      = HERE.parent / "rounds"
FINAL_UNIQUE    = HERE.parent / "final_ranked_unique.csv"
FINAL_RANKED    = HERE.parent / "final_ranked.csv"

WT_PDB_ID       = "1M17"
WT_DIR          = HERE.parent / "_wt_receptor"
SEL_DIR         = HERE.parent / "_selectivity"
SEL_CSV         = HERE.parent / "selectivity.csv"
SEL_HTML        = HERE.parent / "selectivity_view.html"


# --- Helpers ---------------------------------------------------------------

def load_top_unique(csv_path: Path, n: int) -> list[dict]:
    """Read the deduplicated top-N (best score per canonical SMILES)."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        canon_fn = lambda s: (Chem.MolToSmiles(Chem.MolFromSmiles(s), canonical=True)
                              if Chem.MolFromSmiles(s) else s)
    except ImportError:
        canon_fn = lambda s: s

    rows = list(csv.DictReader(csv_path.open()))
    best: dict[str, dict] = {}
    for r in rows:
        c = r.get("canonical_smiles") or canon_fn(r["smiles"])
        if c not in best or float(r["score"]) < float(best[c]["score"]):
            best[c] = r
    return sorted(best.values(), key=lambda r: float(r["score"]))[:n]


def find_ligand_pdbqt(name: str) -> Path | None:
    """Locate the prepared ligand PDBQT — same as docking input."""
    for m in sorted(ROUNDS_DIR.glob(f"round_*/ligprep/{name}.pdbqt"), reverse=True):
        return m
    return None


def find_mut_poses(name: str) -> Path | None:
    for m in sorted(ROUNDS_DIR.glob(f"round_*/docking/poses/{name}_poses.pdbqt"),
                    reverse=True):
        return m
    return None


_MODEL1_RX = re.compile(r"^MODEL\s+1\s*\n(.*?)^ENDMDL",
                        flags=re.MULTILINE | re.DOTALL)

def pdbqt_first_pose_to_pdb(pdbqt_path: Path, out_path: Path) -> Path:
    """Convert MODEL 1 of a Vina poses pdbqt → clean single-model PDB."""
    text = pdbqt_path.read_text(errors="replace")
    m = _MODEL1_RX.search(text)
    body = m.group(1) if m else text
    out_lines = []
    for ln in body.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            row = ln.ljust(80)
            row = "HETATM" + row[6:17] + "LIG" + row[20:]
            out_lines.append(row[:66].rstrip())
    out_lines.append("END")
    out_path.write_text("\n".join(out_lines) + "\n")
    return out_path


def b64_file(p: Path) -> str:
    return base64.b64encode(p.read_bytes()).decode("ascii")


# --- Side-by-side HTML viewer ---------------------------------------------

DUAL_TEMPLATE = """<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>EGFR Selectivity — Mutant vs Wild-type</title>
<style>
  html, body { margin:0; padding:0; height:100%; font-family:-apple-system, sans-serif; background:#0b0c10; color:#e8e8e8; overflow:hidden; }
  header { padding:10px 18px; background:#1a1d23; border-bottom:1px solid #2a2d33; }
  header h1 { margin:0; font-size:16px; color:#6cc7ff; }
  #grid { display:flex; height:calc(100vh - 110px); }
  .pane { flex:1; position:relative; border-left:1px solid #2a2d33; }
  .pane:first-child { border-left:none; }
  .pane h2 { position:absolute; top:8px; left:14px; margin:0; font-size:13px; color:#6cc7ff; z-index:10; background:rgba(0,0,0,0.6); padding:4px 8px; border-radius:4px; }
  .pane .legend { position:absolute; top:36px; left:14px; font-size:11px; color:#ddd; z-index:10; background:rgba(0,0,0,0.6); padding:6px 10px; border-radius:4px; max-width:280px; }
  .row { display:flex; align-items:center; margin:2px 0; }
  .sw { width:10px; height:10px; border-radius:50%; margin-right:6px; display:inline-block; }
  #controls { padding:10px 18px; background:#1a1d23; border-top:1px solid #2a2d33; font-size:13px; display:flex; align-items:center; gap:16px; }
  button { background:#2a2d33; color:#e8e8e8; border:1px solid #3a3d43; padding:5px 12px; border-radius:4px; cursor:pointer; }
  button:hover { background:#3a3d43; }
  select { background:#2a2d33; color:#e8e8e8; border:1px solid #3a3d43; padding:4px 8px; border-radius:4px; }
  .summary { font-size:12px; color:#aaa; }
  .summary b { color:#e8e8e8; }
</style>
</head>
<body>
<header><h1>EGFR Selectivity Check — Mutant (T790M/L858R) vs Wild-type</h1></header>
<div id="grid">
  <div class="pane">
    <h2>Mutant — T790M/L858R (4ZAU)</h2>
    <div class="legend" id="legL"></div>
    <div id="viewL" style="width:100%; height:100%"></div>
  </div>
  <div class="pane">
    <h2>Wild-type (1M17)</h2>
    <div class="legend" id="legR"></div>
    <div id="viewR" style="width:100%; height:100%"></div>
  </div>
</div>
<div id="controls">
  <label>Ligand:
    <select id="ligPick" onchange="showLigand(this.value)"></select>
  </label>
  <button onclick="applyStyle('default')">Default</button>
  <button onclick="applyStyle('surface')">Surface + Stick</button>
  <button onclick="resetViews()">Reset</button>
  <div class="summary" id="summary"></div>
</div>

<!-- Embedded data tags below (receptors + per-ligand pose pairs) -->
__DATA_TAGS__

<script>__THREEDMOL_JS__</script>

<script>
var L = null, R = null;          // viewer handles
var DATA = { ligands: [] };      // populated below

function b64decodeUtf8(b64) {
  var bin = atob(b64);
  var bytes = new Uint8Array(bin.length);
  for (var i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return new TextDecoder('utf-8').decode(bytes);
}

function readReceptor(id) {
  var t = document.getElementById(id);
  return b64decodeUtf8(t.textContent.trim());
}

function buildLigandIndex() {
  var tags = document.querySelectorAll('script[type="chemical/x-base64"][data-role="ligand"]');
  var ligs = {};
  for (var i = 0; i < tags.length; i++) {
    var t = tags[i];
    var name = t.getAttribute('data-name');
    var side = t.getAttribute('data-side');   // 'mut' or 'wt'
    var meta = JSON.parse(t.getAttribute('data-meta') || '{}');
    if (!ligs[name]) ligs[name] = {name: name, meta: meta};
    ligs[name][side] = b64decodeUtf8(t.textContent.trim());
  }
  DATA.ligands = Object.keys(ligs).map(function(k){ return ligs[k]; });
}

function buildPicker() {
  var sel = document.getElementById('ligPick');
  for (var i = 0; i < DATA.ligands.length; i++) {
    var l = DATA.ligands[i];
    var opt = document.createElement('option');
    opt.value = l.name;
    var dd = l.meta.ddG != null ? '  ΔΔG=' + l.meta.ddG.toFixed(2) : '';
    opt.text = l.name + dd;
    sel.appendChild(opt);
  }
}

function styleProtein(v) { v.setStyle({hetflag:false}, {cartoon:{color:'spectrum'}}); }

function paintLigand(v, color) {
  v.setStyle({hetflag:true}, {stick:{colorscheme: color, radius: 0.22}});
}

function showLigand(name) {
  var lig = DATA.ligands.find(function(l){ return l.name === name; });
  if (!lig) return;

  // Mutant pane
  L.removeAllModels();
  L.addModel(readReceptor('mutReceptor'), 'pdb');
  L.addModel(lig.mut, 'pdb');
  styleProtein(L);
  paintLigand(L, 'cyanCarbon');
  L.zoomTo({hetflag:true});
  L.render();

  // Wild-type pane
  R.removeAllModels();
  R.addModel(readReceptor('wtReceptor'), 'pdb');
  R.addModel(lig.wt, 'pdb');
  styleProtein(R);
  paintLigand(R, 'magentaCarbon');
  R.zoomTo({hetflag:true});
  R.render();

  // Per-ligand legend + summary
  var m = lig.meta;
  document.getElementById('legL').innerHTML =
    '<div class="row"><span class="sw" style="background:cyan"></span>' + name +
    '<br>mut score: <b>' + m.mut_score + '</b> kcal/mol</div>';
  document.getElementById('legR').innerHTML =
    '<div class="row"><span class="sw" style="background:magenta"></span>' + name +
    '<br>wt score: <b>' + (m.wt_score != null ? m.wt_score : 'n/a') + '</b> kcal/mol</div>';
  var sel = m.selective === true ? '✓ T790M-selective'
          : m.selective === false ? '✗ non-selective'
          : '— inconclusive';
  document.getElementById('summary').innerHTML =
    '<b>ΔΔG (wt − mut):</b> ' + (m.ddG != null ? m.ddG.toFixed(2) : 'n/a') +
    ' kcal/mol &nbsp;·&nbsp; <b>' + sel + '</b>';
}

var STYLE_MODE = 'default';
function applyStyle(s) {
  STYLE_MODE = s;
  [L, R].forEach(function(v) {
    v.removeAllSurfaces();
    v.setStyle({}, {});
    styleProtein(v);
    paintLigand(v, v === L ? 'cyanCarbon' : 'magentaCarbon');
    if (s === 'surface') {
      v.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.55, color:'white'}, {hetflag:false});
    }
    v.render();
  });
}
function resetViews() {
  [L, R].forEach(function(v) { v.zoomTo({hetflag:true}); v.render(); });
}

window.addEventListener('load', function(){
  if (typeof $3Dmol === 'undefined') {
    document.body.insertAdjacentHTML('beforeend',
      '<div style="padding:24px;color:#ff8a8a">3Dmol.js failed to load.</div>');
    return;
  }
  L = $3Dmol.createViewer(document.getElementById('viewL'), {backgroundColor:'#0b0c10'});
  R = $3Dmol.createViewer(document.getElementById('viewR'), {backgroundColor:'#0b0c10'});
  buildLigandIndex();
  buildPicker();
  if (DATA.ligands.length) showLigand(DATA.ligands[0].name);
});
</script>
</body></html>
"""


# --- Pipeline --------------------------------------------------------------

def prepare_wt() -> dict:
    """Download + clean wild-type EGFR (1M17). Reuse if already prepared."""
    WT_DIR.mkdir(parents=True, exist_ok=True)
    pdbqt = WT_DIR / "receptor.pdbqt"
    fixed = WT_DIR / "protein_prepared.pdb"
    ref_glob = list(WT_DIR.glob("ref_ligand_*.pdb"))
    if pdbqt.exists() and fixed.exists() and ref_glob:
        print(f"  ✓ Wild-type already prepared at {WT_DIR}")
        return {"pdb": fixed, "pdbqt": pdbqt, "ligand": ref_glob[0]}
    print(f"  Preparing wild-type EGFR (PDB {WT_PDB_ID}) …")
    return protprep.prepare_protein(
        source=WT_PDB_ID,
        output_dir=WT_DIR,
        ph=7.4,
        keep_water=False,
        keep_hetero=True,        # keep co-crystal until extraction
        minimize=False,
        extract_ligand_resname="auto",
    )


def dock_into_wt(ligand_pdbqts: list[Path], rec: dict) -> Path:
    SEL_DIR.mkdir(parents=True, exist_ok=True)
    center, size = docking.grid_from_ligand(rec["ligand"], padding=8.0)
    print(f"  WT grid center=({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}) "
          f"size=({size[0]:.0f}, {size[1]:.0f}, {size[2]:.0f}) Å")
    return docking.dock_batch(
        receptor_pdbqt=rec["pdbqt"],
        ligand_files=ligand_pdbqts,
        output_dir=SEL_DIR,
        center=center, size=size,
        backend="vina",
        exhaustiveness=8,
        n_poses=5,
    )


def build_selectivity_csv(top: list[dict], wt_csv: Path) -> list[dict]:
    """Merge mutant scores (from final_ranked_unique.csv) with WT scores."""
    wt_scores = {}
    for r in csv.DictReader(wt_csv.open()):
        try:
            wt_scores[r["ligand"]] = float(r["best_score"])
        except (ValueError, TypeError):
            pass

    out = []
    for r in top:
        name = r["name"]
        mut = float(r["score"])
        wt  = wt_scores.get(name)
        ddG = (wt - mut) if wt is not None else None
        selective = (ddG is not None and ddG >= 1.5)
        out.append({
            "name": name,
            "smiles": r["smiles"],
            "mut_score": round(mut, 3),
            "wt_score":  None if wt is None else round(wt, 3),
            "ddG_wt_minus_mut": None if ddG is None else round(ddG, 3),
            "selective": selective,
        })
    out.sort(key=lambda r: (r["ddG_wt_minus_mut"] is None, -(r["ddG_wt_minus_mut"] or 0)))
    with SEL_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    return out


def build_dual_html(sel_rows: list[dict],
                    mut_rec: Path, wt_rec: Path,
                    mut_poses_dir: Path, wt_poses_dir: Path) -> Path:
    """Build the side-by-side HTML viewer with all ligands embedded."""
    tags: list[str] = []
    tags.append(f'<script type="chemical/x-base64" data-role="receptor" '
                f'id="mutReceptor">{b64_file(mut_rec)}</script>')
    tags.append(f'<script type="chemical/x-base64" data-role="receptor" '
                f'id="wtReceptor">{b64_file(wt_rec)}</script>')

    import json
    for row in sel_rows:
        name = row["name"]
        meta = {
            "mut_score": row["mut_score"],
            "wt_score": row["wt_score"],
            "ddG": row["ddG_wt_minus_mut"],
            "selective": row["selective"],
        }
        mut_pdb = mut_poses_dir / f"{name}_best.pdb"
        wt_pdb  = wt_poses_dir  / f"{name}_best.pdb"
        if not mut_pdb.exists() or not wt_pdb.exists():
            continue
        meta_json = json.dumps(meta).replace('"', "&quot;")
        tags.append(f'<script type="chemical/x-base64" data-role="ligand" '
                    f'data-side="mut" data-name="{name}" '
                    f'data-meta="{meta_json}">{b64_file(mut_pdb)}</script>')
        tags.append(f'<script type="chemical/x-base64" data-role="ligand" '
                    f'data-side="wt" data-name="{name}" '
                    f'data-meta="{meta_json}">{b64_file(wt_pdb)}</script>')

    js = _ensure_3dmol_js() or ""
    html = (DUAL_TEMPLATE
            .replace("__DATA_TAGS__", "\n".join(tags))
            .replace("__THREEDMOL_JS__", js))
    SEL_HTML.write_text(html, encoding="utf-8")
    return SEL_HTML


# --- Main ------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--no-open", action="store_true")
    args = ap.parse_args()

    src_csv = FINAL_UNIQUE if FINAL_UNIQUE.exists() else FINAL_RANKED
    print(f"== Selectivity check ==")
    print(f"Ranked source: {src_csv.name}")

    top = load_top_unique(src_csv, args.top)
    print(f"Top {len(top)} unique hits (best score per canonical SMILES):")

    # Resolve ligprep PDBQTs and mutant poses for cached side-by-side view
    ligand_files: list[Path] = []
    mut_poses_dir = HERE.parent / "_top_pose_cache"
    mut_poses_dir.mkdir(exist_ok=True)
    for r in top:
        name = r["name"]
        lig = find_ligand_pdbqt(name)
        mut_poses = find_mut_poses(name)
        if lig is None or mut_poses is None:
            print(f"  · {name}  ! missing ligprep or mutant pose; skipping")
            continue
        # cache best mutant pose as PDB
        pdbqt_first_pose_to_pdb(mut_poses, mut_poses_dir / f"{name}_best.pdb")
        ligand_files.append(lig)
        print(f"  · {name}  mut_score={r['score']}")

    if not ligand_files:
        sys.exit("No ligands resolvable for selectivity check.")

    print("\n[1/3] Prepare wild-type EGFR (1M17)")
    wt = prepare_wt()
    print(f"  PDB:   {wt['pdb']}")
    print(f"  PDBQT: {wt['pdbqt']}")
    print(f"  Ref:   {wt['ligand']}")

    print(f"\n[2/3] Dock {len(ligand_files)} ligands into wild-type")
    wt_csv = dock_into_wt(ligand_files, wt)
    print(f"  WT docking CSV: {wt_csv}")

    # Cache best WT poses as PDB too, for the dual viewer
    wt_poses_dir = HERE.parent / "_wt_pose_cache"
    wt_poses_dir.mkdir(exist_ok=True)
    for lig in ligand_files:
        name = lig.stem
        wt_pose = SEL_DIR / "poses" / f"{name}_poses.pdbqt"
        if wt_pose.exists():
            pdbqt_first_pose_to_pdb(wt_pose, wt_poses_dir / f"{name}_best.pdb")

    print("\n[3/3] Build selectivity table + comparison viewer")
    sel_rows = build_selectivity_csv(top, wt_csv)

    print(f"\n{'Name':<13} {'mut':>8} {'wt':>8} {'ΔΔG':>8}  selective?")
    print("-" * 50)
    for r in sel_rows:
        mut_s = f"{r['mut_score']:.2f}"
        wt_s  = "n/a " if r["wt_score"] is None else f"{r['wt_score']:.2f}"
        dd_s  = "n/a " if r["ddG_wt_minus_mut"] is None else f"{r['ddG_wt_minus_mut']:+.2f}"
        flag  = "✓" if r["selective"] else "·"
        print(f"{r['name']:<13} {mut_s:>8} {wt_s:>8} {dd_s:>8}  {flag}")

    html = build_dual_html(sel_rows, MUT_RECEPTOR, wt["pdb"], mut_poses_dir, wt_poses_dir)
    print(f"\n✓ selectivity.csv  → {SEL_CSV}")
    print(f"✓ comparison view  → {html}")

    if not args.no_open:
        import subprocess, platform
        if platform.system() == "Darwin":
            subprocess.run(["open", "-a", "Google Chrome", str(html.resolve())])
            print("  Opened in Google Chrome.")


if __name__ == "__main__":
    main()
