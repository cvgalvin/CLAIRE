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


def _split(proj: Project, d: Path) -> None:
    py = tools.python(proj)
    subprocess.run([py, tools.legacy("consolidate_filtered_contacts.py")], cwd=d, check=True)
    out = d / "residue_contact_fuzzballs"
    out.mkdir(exist_ok=True)
    made = {}
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
    proj.write_manifest("contacts", {"residue_fuzzballs": {k: f"residue_contact_fuzzballs/{proj.lig}_{k}_residue_contacts"
                                                           for k in made},
                                     "files": made, "n_outputs": sum(len(v) for v in made.values())})
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
