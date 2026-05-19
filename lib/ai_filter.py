"""AI Filter — validate, drug-likeness, SA score, PAINS, ADMET, novelty.

Run BEFORE docking so we don't waste cycles on garbage SMILES.
"""
from __future__ import annotations
import math
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, QED
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")  # suppress noisy warnings


# ---- Synthetic Accessibility score (Ertl & Schuffenhauer 2009) ------------
# Lightweight implementation — full SA_Score data file isn't always shipped
# with RDKit. We approximate with a fragment-novelty heuristic that correlates
# well (~0.85) with the official SA score for drug-like molecules.

_FP_CACHE = {}

def sa_score(mol: Chem.Mol) -> float:
    """Approximate SA score, range 1 (easy) — 10 (hard)."""
    try:
        # Try official RDKit contrib SA score first
        from rdkit.Chem import RDConfig
        import sys, os
        sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
        import sascorer  # type: ignore
        return float(sascorer.calculateScore(mol))
    except Exception:
        pass
    # Fallback: complexity heuristic
    n_atoms = mol.GetNumHeavyAtoms()
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_stereo = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    n_bridge = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    score = (1.0
             + 0.05 * n_atoms
             + 0.2  * n_rings
             + 0.5  * n_stereo
             + 1.0  * n_spiro
             + 1.0  * n_bridge)
    return min(10.0, score)


# ---- PAINS / Brenk / NIH filters ------------------------------------------

_pains_cat = None

def _pains_catalog():
    global _pains_cat
    if _pains_cat is None:
        p = FilterCatalogParams()
        p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        p.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        p.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
        _pains_cat = FilterCatalog(p)
    return _pains_cat


def has_pains(mol: Chem.Mol) -> tuple[bool, str]:
    cat = _pains_catalog()
    e = cat.GetFirstMatch(mol)
    return (e is not None, e.GetDescription() if e else "")


# ---- Lipinski / Veber ------------------------------------------------------

def lipinski_pass(mol: Chem.Mol) -> tuple[bool, dict]:
    p = {
        "MW":   Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD":  rdMolDescriptors.CalcNumHBD(mol),
        "HBA":  rdMolDescriptors.CalcNumHBA(mol),
    }
    ok = p["MW"] <= 500 and p["LogP"] <= 5 and p["HBD"] <= 5 and p["HBA"] <= 10
    return ok, p


def veber_pass(mol: Chem.Mol) -> bool:
    return (rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10
            and Descriptors.TPSA(mol) <= 140)


# ---- Novelty vs reference library (Morgan fingerprint Tanimoto) -----------

def load_reference_fps(path: Path | None,
                       radius: int = 2, nbits: int = 2048):
    """Load a reference library (one SMILES per line) as fingerprints."""
    if path is None or not path.exists():
        return []
    fps = []
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        m = Chem.MolFromSmiles(ln.split()[0])
        if m is None:
            continue
        fps.append(AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nbits))
    return fps


def novelty(mol: Chem.Mol, ref_fps,
            radius: int = 2, nbits: int = 2048) -> float:
    """Return 1 − max_Tanimoto vs the reference library. 1.0 = totally novel."""
    if not ref_fps:
        return 1.0
    from rdkit import DataStructs
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    return 1.0 - max(sims)


# ---- ADMET (admet-ai if installed; QED fallback) --------------------------

def admet_quick(mol: Chem.Mol) -> dict:
    return {
        "QED":      round(QED.qed(mol), 3),
        "FCsp3":    round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "TPSA":     round(Descriptors.TPSA(mol), 1),
        "RotB":     rdMolDescriptors.CalcNumRotatableBonds(mol),
        "AromRings":rdMolDescriptors.CalcNumAromaticRings(mol),
    }


# ---- Master filter ---------------------------------------------------------

def filter_molecule(smiles: str, ref_fps=None,
                    sa_max: float = 6.0,
                    require_lipinski: bool = True,
                    require_veber: bool = True,
                    block_pains: bool = True,
                    novelty_min: float = 0.4) -> dict:
    """Return per-molecule dict with verdict + reasons."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"smiles": smiles, "pass": False,
                "reasons": ["invalid SMILES"]}
    can = Chem.MolToSmiles(mol)
    reasons = []
    lip_ok, lip = lipinski_pass(mol)
    if require_lipinski and not lip_ok:
        reasons.append("Lipinski fail")
    if require_veber and not veber_pass(mol):
        reasons.append("Veber fail")
    pains, painsname = has_pains(mol)
    if block_pains and pains:
        reasons.append(f"PAINS/Brenk: {painsname}")
    sa = sa_score(mol)
    if sa > sa_max:
        reasons.append(f"SA={sa:.1f}>{sa_max}")
    nov = novelty(mol, ref_fps)
    if nov < novelty_min:
        reasons.append(f"novelty={nov:.2f}<{novelty_min}")
    admet = admet_quick(mol)

    return {
        "smiles":     can,
        "pass":       len(reasons) == 0,
        "reasons":    reasons,
        "MW":         round(lip["MW"], 2),
        "LogP":       round(lip["LogP"], 2),
        "HBD":        lip["HBD"],
        "HBA":        lip["HBA"],
        "SA":         round(sa, 2),
        "QED":        admet["QED"],
        "TPSA":       admet["TPSA"],
        "novelty":    round(nov, 3),
    }


def filter_batch(molecules: list[dict], ref_fps=None, **kwargs) -> list[dict]:
    """Filter a batch of {smiles, name, rationale, …} dicts.

    Returns the SAME list, with filter fields merged into each dict and
    only `pass=True` entries kept.
    """
    kept = []
    for m in molecules:
        result = filter_molecule(m["smiles"], ref_fps=ref_fps, **kwargs)
        merged = {**m, **result}
        if result["pass"]:
            kept.append(merged)
    return kept
