"""Step 4: spatially cluster each per-residue fuzzball (vendored spatial_cluster.py), then, for fragments in
`hbond_fragments`, keep only clusters whose members genuinely hydrogen bond the ligand (ported from workflow.py).

Output: 04_clusters/<Fragment>/<fuzzball>_clustered/<id>_<RES>_cluster_<population>.pdb
        04_clusters/<Fragment>/hbond_clusts/ (for hbond_fragments)
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from .. import tools
from ..executor import Task, run_tasks
from ..project import Project

RESIDUES = ["LYS", "ARG", "ASP", "GLU", "SER", "THR", "GLN", "ASN", "HIS", "TYR", "TRP", "ILE", "LEU", "MET", "PHE", "VAL", "GLY"]
POLAR = ["SER", "THR", "GLN", "ASN", "TYR", "TRP", "GLY", "HIS", "ARG", "LYS", "GLU", "ASP"]


def _frag_key(dirname: str, lig: str) -> str:
    # "<LIG>_Fragment_2_residue_contacts" -> "Fragment_2"
    return dirname[len(lig) + 1:].replace("_residue_contacts", "")


def _fragment_name_map(proj: Project) -> dict:
    """User fragment names -> Fragment_<n>."""
    return {f["name"]: f"Fragment_{i}" for i, f in enumerate(proj.fragments(), 1)}


def _tasks(proj: Project, out: Path) -> list[Task]:
    m = proj.manifest("contacts")
    tasks = []
    for frag, rel in m["residue_fuzzballs"].items():
        src = proj.stage("contacts", False) / rel
        dst = out / frag
        dst.mkdir(parents=True, exist_ok=True)
        for f in sorted(src.glob("*.pdb")):
            if not any(r in f.name for r in RESIDUES):
                continue
            shutil.copy(f, dst / f.name)
            tasks.append(Task(f"{tools.python(proj)} {tools.legacy('spatial_cluster.py')} {f.name} {proj.ligand_pdb()}",
                              str(dst), f"cluster_{frag}_{f.stem}"))
    return tasks


def hbond_screen(proj: Project, fragdir: Path, min_pop: int) -> list[str]:
    """Port of the workflow.py H-bond check: a cluster is kept if >=10% of its members H-bond the ligand."""
    from pyrosetta import Pose, generate_nonstandard_residue_set, init, pose_from_file, rosetta
    init("-pack_missing_sidechains False -mute all")
    kept = []
    for resdir in sorted(p for p in fragdir.iterdir() if p.is_dir() and p.name.endswith("_clustered")):
        parts = resdir.name.split("_")
        rn = parts[3] if len(parts) > 3 else ""
        if rn not in POLAR:
            continue
        for clust in sorted(resdir.glob("*.pdb")):
            pop = int(clust.stem.split("_")[-1])
            p = Pose()
            generate_nonstandard_residue_set(p, [str(proj.params())])
            pose_from_file(p, str(clust))
            lig = p.total_residue()
            nhb = 0
            for i in range(1, lig):
                cp = Pose()
                cp.append_residue_by_jump(p.residue(i).clone(), 1)
                cp.append_residue_by_jump(p.residue(lig).clone(), 1)
                if cp.residue(1).name()[:3] == cp.residue(2).name()[:3]:
                    continue
                hs = rosetta.core.scoring.hbonds.HBondSet()
                cp.update_residue_neighbors()
                rosetta.core.scoring.hbonds.fill_hbond_set(cp, False, hs)
                if hs.nhbonds() > 0:
                    nhb += 1
                if nhb >= 0.1 * pop:
                    kept.append(str(clust))
                    break
    out = fragdir / "hbond_clusts"
    out.mkdir(exist_ok=True)
    final = []
    for c in kept:
        if int(Path(c).stem.split("_")[-1]) >= min_pop:
            shutil.copy(c, out / Path(c).name)
            final.append(Path(c).name)
    print(f"[cluster] {fragdir.name}: {len(final)} H-bonding clusters with population >= {min_pop}")
    return final


def _collect(proj: Project, out: Path) -> None:
    names = _fragment_name_map(proj)
    hb = {names.get(n, n) for n in proj.config.get("hbond_fragments") or []}
    summary = {}
    for fragdir in sorted(p for p in out.iterdir() if p.is_dir() and p.name.startswith("Fragment_")):
        n = sum(1 for _ in fragdir.glob("*_clustered/*.pdb"))
        summary[fragdir.name] = {"n_clusters": n}
        if fragdir.name in hb:
            summary[fragdir.name]["hbond_clusters"] = hbond_screen(
                proj, fragdir, int(proj.config["motifs"]["min_hbond_cluster_population"]))
    proj.write_manifest("cluster", {"fragments": summary, "n_outputs": sum(v["n_clusters"] for v in summary.values())})
    print("[cluster] done. Next: define polar_atoms in claire.yaml, then `claire motifs`.")


def run(proj: Project, args) -> None:
    out = proj.stage("cluster")
    if not getattr(args, "collect", False):
        if not run_tasks(proj, "cluster", _tasks(proj, out), label="cluster"):
            return
    _collect(proj, out)
