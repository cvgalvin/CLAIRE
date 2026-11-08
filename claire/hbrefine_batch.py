"""`claire hbrefine-dir`: run HBRefine (fixed or original) on every PDB in a directory, without a project."""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from importlib import resources
from pathlib import Path

OUTCOMES = {"_unrefined.pdb": "already_satisfied", "_fixed.pdb": "fixed", "_unsats.pdb": "partially_fixed",
            "_ubp.pdb": "ligand_fixed_buried_polar_remain", "_bothproblems.pdb": "partially_fixed"}


AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
      "THR", "TRP", "TYR", "VAL", "HOH"}


def params_name(params: Path) -> str:
    for line in Path(params).read_text().splitlines():
        if line.startswith("NAME "):
            return line.split()[1]
    return ""


def normalize_ligand_name(text: str, name3: str) -> tuple[str, bool]:
    """Rename a ligand whose residue name differs from the params NAME only in case (e.g. 'prg' vs 'PRG'), in the
    coordinate records and in the HETNAM / REMARK 666 records Rosetta reads back."""
    import re
    out, changed = [], False
    for l in text.splitlines(keepends=True):
        new = l
        if l.startswith(("ATOM", "HETATM")) and len(l) > 20:
            rn = l[17:20].strip()
            if rn.upper() not in AA and rn != name3 and rn.upper() == name3.upper():
                new = l[:17] + f"{name3:>3}" + l[20:]
        elif l.startswith(("HETNAM", "REMARK 666")):
            new = re.sub(r"(?<![A-Za-z0-9])" + re.escape(name3) + r"(?![A-Za-z0-9])", name3, l, flags=re.IGNORECASE)
        changed |= new != l
        out.append(new)
    return "".join(out), changed


def _one(pdb: Path, params: Path, out: Path, original: bool, cap: float, python: str, diagnose: bool = False,
         extra: tuple = ()) -> dict:
    d = out / pdb.stem
    d.mkdir(parents=True, exist_ok=True)
    text, changed = normalize_ligand_name(pdb.read_text(), params_name(params))
    if changed:                        # work on a renamed copy; the input file is left untouched
        pdb = d / f"{pdb.stem}.input.pdb"
        pdb.write_text(text)
    if original:
        cmd = [python, str(resources.files("claire.legacy").joinpath("hbrefine.py")), str(params), str(pdb), str(d)]
    else:
        cmd = [python, "-m", "claire.hbrefine_fixed", str(params), str(pdb), str(d), str(cap)] + (["diagnose"] if diagnose else []) + list(extra)
    with open(d / "hbrefine.log", "w") as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    stem = d.name
    row = {"design": stem, "version": "original" if original else f"fixed (cap {cap} REU)"}
    js = d / f"{pdb.stem}_hbrefine.json"
    if js.exists():
        s = json.loads(js.read_text())
        hs = s.get("hbond_search") or []
        row.update(hbond_candidates=sum(h["candidates"] for h in hs), hbond_reducing=sum(h["reducing"] for h in hs),
                   hbond_reducing_within_cap=sum(h["reducing_within_cap"] for h in hs),
                   best_reducing_dE=min((h["best_reducing_dE"] for h in hs if h["best_reducing_dE"] is not None), default=None),
                   best_reducing_dE_fa_rep=next((h.get("best_reducing_dE_fa_rep") for h in hs if h.get("best_reducing_dE_fa_rep") is not None), None),
                   best_reducing_dE_after_min=min((h["best_reducing_dE_after_min"] for h in hs
                                                   if h.get("best_reducing_dE_after_min") is not None), default=None))
        h0 = next((h for h in hs if h.get("control_self_dE") is not None), {})
        row.update(best_min_mutation=h0.get("best_min_mutation"),
                   **{k: h0.get(k) for k in ("control_self_mutation", "control_self_dE", "control_self_dE_after_min",
                                             "control_np_mutation", "control_np_dE", "control_np_dE_after_min")})
        row.update(outcome=s.get("outcome"), mutations=" ".join(s.get("mutations", [])),
                   n_mutations=len(s.get("mutations", [])), energy_change=s.get("energy_change"),
                   unsat_before=" ".join(s.get("unsat_before", [])), unsat_after=" ".join(s.get("unsat_after", [])))
    else:
        outs = sorted(p.name for p in d.glob("*.pdb"))
        outcome = next((v for k, v in OUTCOMES.items() for o in outs if o.endswith(k)), "no_solution" if rc == 0 else "error")
        nm = next((o.split("mutations")[0].rsplit("_", 1)[-1] for o in outs if "mutations" in o), "0")
        row.update(outcome=outcome, n_mutations=nm)
    if rc != 0 and row.get("outcome") in (None, "no_solution"):
        row["outcome"] = "error"
    row["output"] = " ".join(p.name for p in d.glob("*.pdb") if not p.name.endswith(".input.pdb"))
    return row


def run(params: str, indir: str, outdir: str, original: bool, cap: float, workers: int, python: str | None,
        diagnose: bool = False, extra: tuple = ()) -> int:
    params, indir, outdir = Path(params).resolve(), Path(indir).resolve(), Path(outdir).resolve()
    pdbs = sorted(indir.glob("*.pdb"))
    if not pdbs:
        print(f"No .pdb files in {indir}", file=sys.stderr)
        return 1
    outdir.mkdir(parents=True, exist_ok=True)
    python = python or sys.executable
    print(f"HBRefine ({'original' if original else f'fixed, cap {cap} REU'}) on {len(pdbs)} structures, {workers} workers")
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_one, p, params, outdir, original, cap, python, diagnose, extra): p for p in pdbs}
        for i, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if i % max(1, len(pdbs) // 20) == 0 or i == len(pdbs):
                print(f"  {i}/{len(pdbs)}", flush=True)
    rows.sort(key=lambda r: r["design"])
    keys = ["design", "version", "outcome", "n_mutations", "mutations", "energy_change", "unsat_before", "unsat_after",
            "hbond_candidates", "hbond_reducing", "hbond_reducing_within_cap", "best_reducing_dE", "best_reducing_dE_fa_rep",
            "best_reducing_dE_after_min", "best_min_mutation", "control_self_mutation", "control_self_dE",
            "control_self_dE_after_min", "control_np_mutation", "control_np_dE", "control_np_dE_after_min", "output"]
    with open(outdir / "hbrefine_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    counts = {}
    for r in rows:
        counts[r.get("outcome")] = counts.get(r.get("outcome"), 0) + 1
    print(f"Outcomes: {counts}\nSummary: {outdir / 'hbrefine_summary.csv'}; per-design outputs and logs in {outdir}/<design>/")
    return 0
