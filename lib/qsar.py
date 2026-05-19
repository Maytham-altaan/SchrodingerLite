"""Strike — QSAR / activity modeling (Schrödinger Strike equivalent)."""
from __future__ import annotations
import argparse
import csv
import pickle
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .utils import banner, new_project


def _fp_matrix(smiles: list[str], radius: int = 2, nbits: int = 2048):
    X = np.zeros((len(smiles), nbits), dtype=np.uint8)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nbits)
        arr = np.zeros(nbits, dtype=np.uint8)
        from rdkit import DataStructs
        DataStructs.ConvertToNumpyArray(fp, arr)
        X[i] = arr
    return X


def train(csv_path: Path, smiles_col: str, y_col: str,
          output_dir: Path, task: str = "regression",
          model_type: str = "rf") -> dict:
    from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
    from sklearn.model_selection import cross_val_score, train_test_split
    from sklearn.metrics import r2_score, mean_absolute_error, roc_auc_score

    output_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(csv_path.open()))
    smiles = [r[smiles_col] for r in rows]
    y = np.array([float(r[y_col]) for r in rows])
    X = _fp_matrix(smiles)

    if task == "regression":
        if model_type == "xgb":
            from xgboost import XGBRegressor
            model = XGBRegressor(n_estimators=400, max_depth=6, n_jobs=-1)
        else:
            model = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=42)
    else:
        if model_type == "xgb":
            from xgboost import XGBClassifier
            model = XGBClassifier(n_estimators=400, max_depth=6, n_jobs=-1)
        else:
            model = RandomForestClassifier(n_estimators=500, n_jobs=-1, random_state=42)

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=42)
    model.fit(Xtr, ytr)
    pred = model.predict(Xte)

    metrics = {}
    if task == "regression":
        metrics["R2"]  = round(r2_score(yte, pred), 3)
        metrics["MAE"] = round(mean_absolute_error(yte, pred), 3)
    else:
        try:
            proba = model.predict_proba(Xte)[:,1]
            metrics["AUC"] = round(roc_auc_score(yte, proba), 3)
        except Exception:
            metrics["accuracy"] = round((pred == yte).mean(), 3)

    model_path = output_dir / "model.pkl"
    with open(model_path, "wb") as f:
        pickle.dump({"model": model, "task": task,
                     "smiles_col": smiles_col, "y_col": y_col}, f)
    return {"metrics": metrics, "model": model_path,
            "n_train": len(Xtr), "n_test": len(Xte)}


def predict(model_path: Path, smiles_list: list[str]) -> list[float]:
    with open(model_path, "rb") as f:
        bundle = pickle.load(f)
    model = bundle["model"]
    X = _fp_matrix(smiles_list)
    return model.predict(X).tolist()


def main(argv=None):
    p = argparse.ArgumentParser(prog="strike",
        description="Strike — QSAR / activity prediction.")
    sub = p.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("train", help="Train a model from CSV")
    t.add_argument("csv", type=Path, help="CSV with SMILES + activity columns")
    t.add_argument("--smiles-col", default="SMILES")
    t.add_argument("--y-col", required=True)
    t.add_argument("--task", choices=["regression","classification"], default="regression")
    t.add_argument("--model", choices=["rf","xgb"], default="rf")
    t.add_argument("-o","--output", type=Path, default=None)

    pr = sub.add_parser("predict", help="Predict activities for new SMILES")
    pr.add_argument("model", type=Path)
    pr.add_argument("smiles", nargs="+")

    args = p.parse_args(argv)

    if args.cmd == "train":
        banner("Strike — QSAR training", f"{args.task} ({args.model})")
        out = new_project("strike") / "output" if args.output is None else args.output
        res = train(args.csv, args.smiles_col, args.y_col, out,
                    task=args.task, model_type=args.model)
        print(f"\n✓ Trained on {res['n_train']} mols, tested on {res['n_test']}")
        print(f"  Metrics: {res['metrics']}")
        print(f"  Model saved: {res['model']}")
    else:
        banner("Strike — prediction")
        preds = predict(args.model, args.smiles)
        for s, p in zip(args.smiles, preds):
            print(f"  {s}\t{p:.3f}")


if __name__ == "__main__":
    main()
