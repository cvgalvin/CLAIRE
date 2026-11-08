"""Step 3: scrape protein contacts with each fragment from PDB structures, filter, and split by residue type.

Phases (run in order; on a cluster each submitted phase is resumed with `claire contacts --collect`):
  extract  -> download structures, align contacts onto each fragment (vendored extract_align_contacts.py)
  filter   -> redundancy/quality filtering in chunks (vendored filter_contacts_array.py)
  split    -> consolidate fuzzballs and split into per-residue fuzzballs (vendored scripts)
Output: 03_contacts/residue_contact_fuzzballs/<LIG>_<Fragment>_residue_contacts/*.pdb
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

from .. import sources, tools
from ..executor import Task, run_tasks
from ..project import ClaireError, Project

CHUNK_PDBS = 150       # structures per extraction task (as in process_fragment_search.py)
CHUNK_FILTER = 500     # contact files per filtering task (as in setup_filter_contacts.py)


def _state(d: Path) -> dict:
    p = d / "state.json"
    return json.loads(p.read_text()) if p.exists() else {"phase": "extract"}


def _save(d: Path, st: dict) -> None:
    (d / "state.json").write_text(json.dumps(st, indent=1))


def _extract_tasks(proj: Project, d: Path) -> list[Task]:
    res = json.loads((proj.stage("fragments", False) / "PDB_search_results.json").read_text())
    cap = proj.config["sources"].get("max_entries_per_fragment")
    if cap:
        import random
        for frag, v in res.items():
            if len(v["PDBs"]) > int(cap):
                n0 = len(v["PDBs"])
                order = sorted(v["PDBs"])
                random.Random(int(proj.config["sources"].get("subsample_seed", 1))).shuffle(order)
                v["PDBs"] = sorted(order[:int(cap)])     # nested: the first 1000 are also in a 2000-entry subset
                print(f"[contacts] {frag}: using a random subset of {cap} of {n0} PDB entries (sources.max_entries_per_fragment)")
        (d / "entries_used.json").write_text(json.dumps({k: v["PDBs"] for k, v in res.items()}))
    all_ids = sorted({p for v in res.values() for p in v["PDBs"]})
    print(f"[contacts] fetching {len(all_ids)} structures "
          f"({proj.config['sources']['structures']}, fallback to RCSB: {proj.config['sources']['fallback_to_rcsb']})")
    root, got = sources.fetch_structures(proj, all_ids)
    counts = {}
    for v in got.values():
        counts[v] = counts.get(v, 0) + 1
    print(f"[contacts] structures available in {root}: {counts}")
    (d / "structure_sources.json").write_text(json.dumps(got, indent=0))
    inp = d / "Alignment_Inputs"
    shutil.rmtree(inp, ignore_errors=True)
    inp.mkdir()
    tasks, n = [], 0
    env = (f"CLAIRE_PDB_DIR={root} CLAIRE_FRAGMENT_DIR={proj.stage('fragments', False)} ")
    for frag, v in res.items():
        ids = [p for p in v["PDBs"] if got.get(p.upper()) != "missing"]
        for i in range(0, len(ids), CHUNK_PDBS):
            n += 1
            f = inp / f"PDB_search_results_{n}.json"
            f.write_text(json.dumps({frag: {"PDBs": ids[i:i + CHUNK_PDBS], "Ligands": v["Ligands"], "SMILES": v["SMILES"]}}))
            tasks.append(Task(f"{env}{tools.python(proj)} {tools.legacy('extract_align_contacts.py')} {f}",
                              str(d), f"extract_{n}"))
    return tasks


def _filter_tasks(proj: Project, d: Path) -> list[Task]:
    tap = d / "Transformed_Aligned_PDBs"
    if not tap.exists():
        raise ClaireError("No aligned contacts were produced; check 03_contacts/logs/extract_*.log")
    dirs = []
    for fd in sorted(p for p in tap.iterdir() if p.is_dir()):
        files = sorted(f for f in fd.iterdir() if f.suffix == ".pdb")
        if len(files) <= CHUNK_FILTER:
            dirs.append(fd)
        else:
            for k, i in enumerate(range(0, len(files), CHUNK_FILTER), 1):
                sub = fd / "filter_input" / str(k)
                sub.mkdir(parents=True, exist_ok=True)
                for f in files[i:i + CHUNK_FILTER]:
                    shutil.copy2(f, sub / f.name)
                dirs.append(sub)
    return [Task(f"{tools.python(proj)} {tools.legacy('filter_contacts_array.py')} {p}", str(d), f"filter_{i}")
            for i, p in enumerate(dirs, 1)]


def drop_impossible_contacts(path: Path, fragment_atoms: set, lig: str, cutoff_polar: float = 2.4,
                             cutoff_other: float = 2.8) -> tuple[int, int]:
    """Remove contact residues that approach the fragment closer than physically possible.

    A residue is removed if any of its heavy atoms is closer to a fragment heavy atom than `cutoff_polar` (A) for an
    N/O...N/O pair (shorter than all but the strongest hydrogen bonds), or `cutoff_other` (A) for any pair involving C or
    S (well inside van der Waals contact). Such contacts come from extraction or alignment artefacts. Returns
    (residues kept, residues removed); the file is rewritten in place.
    """
    import numpy as np
    lines = path.read_text().splitlines(keepends=True)

    def element(l):
        el = l[76:78].strip().upper() if len(l) > 77 else ""
        return el or l[12:16].strip()[:1].upper()

    def is_lig(l):
        return l[21] == "X" or l[17:20].strip().upper() == lig.upper()

    frag = [(element(l), (float(l[30:38]), float(l[38:46]), float(l[46:54]))) for l in lines
            if l.startswith(("ATOM", "HETATM")) and is_lig(l) and element(l) != "H" and l[12:16].strip() in fragment_atoms]
    if not frag:
        return 0, 0
    fxyz = np.array([x for _, x in frag])
    fpolar = np.array([e in ("N", "O") for e, _ in frag])
    res = {}
    for l in lines:
        if l.startswith(("ATOM", "HETATM")) and not is_lig(l):
            res.setdefault((l[21], l[22:27]), []).append(l)
    bad = set()
    for key, rl in res.items():
        for l in rl:
            e = element(l)
            if e == "H":
                continue
            d = np.linalg.norm(fxyz - np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), axis=1)
            lim = np.where(fpolar & (e in ("N", "O")), cutoff_polar, cutoff_other)
            if (d < lim).any():
                bad.add(key)
                break
    if bad:
        path.write_text("".join(l for l in lines
                                if not (l.startswith(("ATOM", "HETATM", "ANISOU", "TER")) and not is_lig(l)
                                        and (l[21], l[22:27]) in bad)))
    return len(res) - len(bad), len(bad)

def _split(proj: Project, d: Path) -> None:
    py = tools.python(proj)
    subprocess.run([py, tools.legacy("consolidate_filtered_contacts.py")], cwd=d, check=True)
    out = d / "residue_contact_fuzzballs"
    out.mkdir(exist_ok=True)
    made, removed_by = {}, {}
    frags = proj.fragments()
    ccfg = proj.config.get("contacts") or {}
    cut_p, cut_o = float(ccfg.get("min_distance_polar", 2.4)), float(ccfg.get("min_distance_other", 2.8))
    for fd in sorted(p for p in (d / "Transformed_Aligned_PDBs").iterdir() if p.is_dir()):
        fb = fd / f"{fd.name}_contact_fuzzball.pdb"
        if not fb.exists():
            print(f"[contacts] no fuzzball for {fd.name}; skipped")
            continue
        subprocess.run([py, tools.legacy("process_filtered_contacts2.py"), str(fb), str(proj.ligand_pdb()),
                        proj.lig, fd.name], cwd=d, check=True)
        rd = d / f"{proj.lig}_{fd.name}_residue_contacts"
        dst = out / rd.name
        shutil.rmtree(dst, ignore_errors=True)
        shutil.move(str(rd), dst)
        made[fd.name] = sorted(p.name for p in dst.glob("*.pdb"))
        frag_atoms = set(frags[int(fd.name.split("_")[1]) - 1]["atoms"]) if fd.name.startswith("Fragment_") else set()
        kept = removed = 0
        for f in dst.glob("*.pdb"):
            k, r = drop_impossible_contacts(f, frag_atoms, proj.lig, cut_p, cut_o)
            kept, removed = kept + k, removed + r
        removed_by[fd.name] = {"kept": kept, "removed": removed}
        if removed:
            print(f"[contacts] {fd.name}: removed {removed} of {kept + removed} contacts closer to the fragment than "
                  f"{cut_o} A (C/S pairs) or {cut_p} A (N/O pairs): overlapping, not genuine contacts")
    proj.write_manifest("contacts", {"residue_fuzzballs": {k: f"residue_contact_fuzzballs/{proj.lig}_{k}_residue_contacts"
                                                           for k in made},
                                     "files": made, "n_outputs": sum(len(v) for v in made.values()),
                                     "overlap_filter": {"min_distance_polar": cut_p, "min_distance_other": cut_o,
                                                        "per_fragment": removed_by}})
    print(f"[contacts] per-residue fuzzballs in {out}. Next: `claire cluster`.")


def run(proj: Project, args) -> None:
    d = proj.stage("contacts")
    st = _state(d) if getattr(args, "collect", False) else {"phase": "extract"}
    while True:
        phase = st["phase"]
        if phase == "extract":
            done = run_tasks(proj, "contacts", _extract_tasks(proj, d), label="extract")
            st["phase"] = "filter"
        elif phase == "filter":
            done = run_tasks(proj, "contacts", _filter_tasks(proj, d), label="filter")
            st["phase"] = "split"
        elif phase == "split":
            _split(proj, d)
            st["phase"] = "done"
            _save(d, st)
            return
        else:
            print("[contacts] already complete.")
            return
        _save(d, st)
        if not done:
            return


def refilter(proj: Project, args) -> None:
    """Apply the overlap filter to an existing project's residue contact files (then rerun `claire cluster`)."""
    d = proj.stage("contacts", False)
    out = d / "residue_contact_fuzzballs"
    if not out.exists():
        raise ClaireError("No residue contact files; run `claire contacts` first.")
    ccfg = proj.config.get("contacts") or {}
    cut_p = float(getattr(args, "cutoff_polar", None) or ccfg.get("min_distance_polar", 2.4))
    cut_o = float(getattr(args, "cutoff_other", None) or ccfg.get("min_distance_other", 2.8))
    frags = proj.fragments()
    summary = {}
    for rd in sorted(out.glob(f"{proj.lig}_Fragment_*_residue_contacts")):
        fname = rd.name[len(proj.lig) + 1:-len("_residue_contacts")]
        atoms = set(frags[int(fname.split("_")[1]) - 1]["atoms"])
        kept = removed = 0
        for f in rd.glob("*.pdb"):
            k, r = drop_impossible_contacts(f, atoms, proj.lig, cut_p, cut_o)
            kept, removed = kept + k, removed + r
        summary[fname] = {"kept": kept, "removed": removed}
        print(f"[contacts] {fname}: removed {removed} of {kept + removed} contacts (C/S pairs < {cut_o} A, N/O pairs < {cut_p} A)")
    m = proj.manifest("contacts")
    m["overlap_filter"] = {"min_distance_polar": cut_p, "min_distance_other": cut_o, "per_fragment": summary}
    proj.write_manifest("contacts", m)
    print("Next: `claire cluster` (re-clusters from the filtered contacts).")
