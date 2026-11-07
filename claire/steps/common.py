"""Helpers shared by workflow steps: score files, threshold filtering, file collection."""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from ..project import ClaireError


def read_scores(paths) -> list[dict]:
    """Read Rosetta JSON score files (one JSON object per line)."""
    scores = []
    for p in paths:
        for line in Path(p).read_text().splitlines():
            line = line.strip()
            if line:
                try:
                    scores.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return scores


ALIASES = {"lighphobesasa": "lig_hphobe_sasa"}   # names used in the original workflow.py -> names in the XMLs


def apply_filters(scores: list[dict], rules: list, label: str = "", skip_missing: bool = False) -> list[dict]:
    """Keep score dicts passing all [metric, op, threshold] rules ("<=" or ">=", as in the original return_filtered)."""
    present = set().union(*[s.keys() for s in scores]) if scores else set()
    kept = scores
    print(f"[filter{(' ' + label) if label else ''}] {len(scores)} designs")
    if not scores:
        return []
    for metric, op, thr in rules:
        name = metric if metric in present else ALIASES.get(metric, metric)
        if name not in present and skip_missing:
            print(f"  {metric}: not computed at this stage, applied later")
            continue
        if name not in present:
            avail = sorted(k for k in present if k not in ("decoy", "nstruct"))
            raise ClaireError(f"Filter metric '{metric}' is not in the score file ({label}). Available metrics: "
                              f"{', '.join(avail)}. Fix the name under filters in claire.yaml.")
        def ok(d, name=name):
            try:
                v = float(d[name])
            except (KeyError, TypeError, ValueError):
                return False
            return v <= float(thr) if op == "<=" else v >= float(thr)
        kept = [d for d in kept if ok(d)]
        print(f"  {name} {op} {thr}: {len(kept)} remain")
    return kept


def decoy_pdb(decoy: str, search_dirs) -> Path | None:
    name = decoy if decoy.endswith(".pdb") else decoy + ".pdb"
    for d in search_dirs:
        p = Path(d) / name
        if p.exists():
            return p
    return None


def link_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.symlink(os.path.abspath(src), dst)
    except OSError:
        shutil.copy(src, dst)


def pdbs(d: Path) -> list[Path]:
    return sorted(p for p in Path(d).glob("*.pdb"))


def require(paths, what: str):
    if not paths:
        raise ClaireError(f"No {what} found.")
    return paths


def match_decoy(decoy: str, stems: dict):
    """Map a Rosetta score-file decoy name to an input file (handles '.pdb' and appended '_0001')."""
    d = decoy[:-4] if decoy.endswith(".pdb") else decoy
    for cand in (d, d[:-5] if d.endswith("_0001") else None):
        if cand and cand in stems:
            return stems[cand]
    return None
