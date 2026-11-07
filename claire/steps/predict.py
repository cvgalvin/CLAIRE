"""Steps 13-14: ColabFold structure prediction (single-sequence mode) and final selection by pLDDT and Ca RMSD.

The analysis ports workflow.py section 15: mean per-residue pLDDT of each predicted model and Ca RMSD between
the prediction and the design model; a design is accepted if any model passes both thresholds, and the model
with the best pLDDT (among passing models) is reported.
"""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np

from .. import tools
from ..executor import Task, run_tasks
from ..project import ClaireError, Project
from .common import link_or_copy

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V"}
BATCH = 50


def _ca(pdb: Path):
    seq, ca, b = [], [], []
    for l in Path(pdb).read_text().splitlines():
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            seq.append(THREE.get(l[17:20], "X"))
            ca.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            b.append(float(l[60:66]))
    return "".join(seq), np.array(ca), np.array(b)


def ca_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    n = min(len(a), len(b))
    a, b = a[:n] - a[:n].mean(0), b[:n] - b[:n].mean(0)
    u, s, vt = np.linalg.svd(a.T @ b)
    dd = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, dd]) @ vt
    return float(np.sqrt(((a @ r - b) ** 2).sum(1).mean()))


def predict(proj: Project, args) -> None:
    d = proj.stage("predict")
    src = proj.stage("fastdesign", False) / "passed"
    designs = [src / n for n in proj.manifest("fastdesign")["designs"]]
    if not getattr(args, "collect", False):
        cb = tools.executable(proj, "colabfold_batch")
        c = proj.config["predict"]
        tasks = []
        for k in range(0, len(designs), BATCH):
            bdir = d / f"batch_{k // BATCH + 1}"
            fdir = bdir / "fastas"
            fdir.mkdir(parents=True, exist_ok=True)
            for p in designs[k:k + BATCH]:
                seq, _, _ = _ca(p)
                (fdir / f"{p.stem}.fasta").write_text(f">{p.stem}\n{seq}\n")
            cmd = (f"{cb} --msa-mode {c['msa_mode']} --num-models {c['num_models']} --num-recycle {c['num_recycles']} "
                   f"--model-order {c['model_order']} --stop-at-score 100 {fdir} {bdir / 'output'}")
            tasks.append(Task(cmd, str(bdir), f"colabfold_{k // BATCH + 1}"))
        if not run_tasks(proj, "predict", tasks, label="ColabFold"):
            return
    preds = sorted(d.glob("batch_*/output/*_unrelaxed_*.pdb"))
    proj.write_manifest("predict", {"n_outputs": len(preds)})
    print(f"[predict] {len(preds)} predicted models. Next: `claire final`.")


def final(proj: Project, args) -> None:
    d = proj.stage("final")
    src = proj.stage("fastdesign", False) / "passed"
    pdir = proj.stage("predict", False)
    cfg = proj.config["final"]
    rows, best = [], {}
    for n in proj.manifest("fastdesign")["designs"]:
        stem = Path(n).stem
        _, dca, _ = _ca(src / n)
        for m in sorted(pdir.glob(f"batch_*/output/{stem}_unrelaxed_*.pdb")):
            _, pca, b = _ca(m)
            plddt = float(b.mean()) if len(b) else float("nan")
            rmsd = ca_rmsd(dca, pca)
            ok = plddt >= float(cfg["plddt_min"]) and rmsd <= float(cfg["ca_rmsd_max"])
            rows.append({"design": stem, "model": m.name, "plddt": round(plddt, 2), "ca_rmsd": round(rmsd, 3), "pass": ok})
            if ok and (stem not in best or plddt > best[stem]["plddt"]):
                best[stem] = {"model": m, "plddt": plddt, "ca_rmsd": rmsd}
    if not rows:
        raise ClaireError("No ColabFold models found; run `claire predict` (and `--collect` on a cluster).")
    with open(d / "all_models.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    out = d / "designs"
    out.mkdir(exist_ok=True)
    with open(d / "final_designs.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["design", "sequence", "best_model", "plddt", "ca_rmsd"])
        for stem, v in sorted(best.items(), key=lambda kv: -kv[1]["plddt"]):
            link_or_copy(src / f"{stem}.pdb", out / f"{stem}.pdb")
            shutil.copy(v["model"], out / f"{stem}_colabfold.pdb")
            seq, _, _ = _ca(src / f"{stem}.pdb")
            w.writerow([stem, seq, v["model"].name, round(v["plddt"], 2), round(v["ca_rmsd"], 3)])
    proj.write_manifest("final", {"designs": sorted(best), "n_outputs": len(best)})
    print(f"[final] {len(best)} designs pass pLDDT >= {cfg['plddt_min']} and Ca RMSD <= {cfg['ca_rmsd_max']} A. "
          f"See {d / 'final_designs.csv'}.")
