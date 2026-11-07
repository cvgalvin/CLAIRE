"""Step 5: assemble binding motifs (Rosetta match constraint files).

1. Polar motifs: one idealized H-bond constraint per ligand polar atom (config `polar_atoms`), with partner
   residues sampled under the original rules (ported verbatim from workflow.py, section 5).
2. Nonpolar contacts: the most populated clusters are scored (fa_atr/fa_rep) and converted to constraint
   blocks from their observed geometry, weighted by cluster population.
3. Hybrid motifs: each polar motif combined with 1 or 2 nonpolar contacts -> 05_motifs/1np, 05_motifs/2np.
"""
from __future__ import annotations

import math
import random
from pathlib import Path

import numpy as np

from ..project import ClaireError, Project

Q_POLAR = ["LYS", "ARG", "SER", "THR", "TYR", "TRP", "GLY"]   # partners for ligand acceptors, charged allowed
POLAR = ["SER", "THR", "TYR", "TRP", "GLY"]
PROBLEMATIC = ["LYS", "ARG", "GLN", "ASN"]
ACCEPTORS = ["GLU", "ASP", "GLY"]                                # partners for ligand donors
HYB = {"sp2": "y", "sp3": "n", "sp": "s"}


# ------------------------------------------------------------------ constraint blocks (from workflow.py)
def polar_block(distanceAB, angleA, angleB, torsionA, torsionAB, torsionB, residue_resname, ligand_resname,
                ligand_atoms, residue_atoms, distance_tolerance=0.15, angle_A_tolerance=15, angle_B_tolerance=15,
                torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                distance_constraint_sample_number=0, isbb=False):
    tag = "residue3"
    head = ["  TEMPLATE::   ATOM_MAP: 1 atom_name: {}".format(" ".join(ligand_atoms)),
            "  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n".format(ligand_resname),
            "  TEMPLATE::   ATOM_MAP: 2 atom_type: {}".format(residue_atoms)]
    if isbb:
        head.append("  TEMPLATE::   ATOM_MAP: 2 is_backbone")
    head.append("  TEMPLATE::   ATOM_MAP: 2 {}: {}\n".format(tag, residue_resname))
    return head + [
        "  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}".format(distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
        "  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}".format(float(angleA), angle_A_tolerance, 100, angle_constraint_sample_number),
        "  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}".format(float(angleB), angle_B_tolerance, 100, angle_constraint_sample_number),
        "  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  180.00  {3:3}".format(float(torsionA), torsion_A_tolerance, 100, torsion_constraint_sample_number),
        "  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}".format(float(torsionB), torsion_B_tolerance, 100, torsion_constraint_sample_number),
        "  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}".format(float(torsionAB), torsion_AB_tolerance, 100, torsion_constraint_sample_number)]


def nonpolar_block(distanceAB, angleA, angleB, torsionA, torsionAB, torsionB, residue_resname, ligand_resname,
                   ligand_atoms, residue_atoms, distance_tolerance=0.25, angle_A_tolerance=15, angle_B_tolerance=15,
                   torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                   torsion_constraint_sample_number=2, angle_constraint_sample_number=2, distance_constraint_sample_number=2):
    return [
        "  TEMPLATE::   ATOM_MAP: 1 atom_name: {}".format(" ".join(ligand_atoms)),
        "  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n".format(ligand_resname),
        "  TEMPLATE::   ATOM_MAP: 2 atom_name: {}".format(" ".join(residue_atoms)),
        "  TEMPLATE::   ATOM_MAP: 2 {}: {}\n".format("residue3", residue_resname),
        "  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}".format(distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
        "  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}".format(float(angleA), angle_A_tolerance, 100, angle_constraint_sample_number),
        "  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}".format(float(angleB), angle_B_tolerance, 100, angle_constraint_sample_number),
        "  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}".format(float(torsionA), torsion_A_tolerance, 100, torsion_constraint_sample_number),
        "  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}".format(float(torsionB), torsion_B_tolerance, 100, torsion_constraint_sample_number),
        "  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}".format(float(torsionAB), torsion_AB_tolerance, 100, torsion_constraint_sample_number)]


# config key -> keyword of polar_block / nonpolar_block
OVERRIDE_KEYS = {
    "distance_tolerance": "distance_tolerance", "angle_A_tolerance": "angle_A_tolerance",
    "angle_B_tolerance": "angle_B_tolerance", "torsion_A_tolerance": "torsion_A_tolerance",
    "torsion_AB_tolerance": "torsion_AB_tolerance", "torsion_B_tolerance": "torsion_B_tolerance",
    "distance_samples": "distance_constraint_sample_number", "angle_samples": "angle_constraint_sample_number",
    "torsion_samples": "torsion_constraint_sample_number",
}
_OVR = {"hbond": {}, "nonpolar": {}}     # set by run() from claire.yaml / --fast


def set_overrides(hbond: dict | None, nonpolar: dict | None) -> None:
    for kind, d in (("hbond", hbond), ("nonpolar", nonpolar)):
        d = d or {}
        bad = [k for k in d if k not in OVERRIDE_KEYS]
        if bad:
            raise ClaireError(f"Unknown motifs.constraints.{kind} key(s) {bad}; allowed: {sorted(OVERRIDE_KEYS)}")
        _OVR[kind] = {OVERRIDE_KEYS[k]: v for k, v in d.items() if v is not None}


def _apply(kind: str, kw: dict) -> dict:
    """Apply user overrides; tolerances of 360 (unconstrained degrees of freedom in workflow.py) stay unconstrained."""
    for k, v in _OVR[kind].items():
        if k.endswith("_tolerance") and kw.get(k) == 360:
            continue
        kw[k] = v
    return kw


def _hbond_block(res, lig, lig_atoms, donor, hyb):
    """Exact parameter sets from workflow.py for each (ligand donor/acceptor, hybridization, Gly/side chain) case,
    then any user overrides (motifs.constraints.hbond / --fast)."""
    if not donor:
        isbb = res == "GLY"
        kw = dict(residue_resname=res, ligand_resname=lig, ligand_atoms=lig_atoms, residue_atoms="HNbb" if isbb else "Hpol",
                  distance_tolerance=0.5, angle_A_tolerance=10, angle_B_tolerance=15, torsion_AB_tolerance=360,
                  torsion_B_tolerance=360, torsion_constraint_sample_number=5, angle_constraint_sample_number=1, isbb=isbb)
        if hyb == "y":
            geo, kw2 = (1.9, 120, 165, 180, 0, 0), dict(torsion_A_tolerance=20, distance_constraint_sample_number=3)
        elif hyb == "s":
            geo, kw2 = (1.9, 180, 165, 180, 0, 0), dict(torsion_A_tolerance=360, distance_constraint_sample_number=3)
        else:
            geo, kw2 = (1.9, 120, 180, 180, 0, 0), dict(torsion_A_tolerance=360, distance_constraint_sample_number=5)
        return polar_block(*geo, **_apply("hbond", {**kw, **kw2}))
    isbb = res == "GLY"
    kw = dict(residue_resname=res, ligand_resname=lig, ligand_atoms=lig_atoms, residue_atoms="OCbb" if isbb else "OOC",
              distance_tolerance=0.5, angle_A_tolerance=15, angle_B_tolerance=10, torsion_A_tolerance=360,
              torsion_AB_tolerance=360, torsion_B_tolerance=20, torsion_constraint_sample_number=5,
              angle_constraint_sample_number=1, distance_constraint_sample_number=3, isbb=isbb)
    return polar_block(1.9, 165, 120, 0, 0, 180, **_apply("hbond", kw))


def polar_motifs(proj: Project) -> list[list[list[str]]]:
    spec = proj.config.get("polar_atoms") or []
    if not spec:
        raise ClaireError("No polar_atoms defined in claire.yaml (one entry per ligand polar atom to H-bond).")
    m = proj.config["motifs"]
    lig = proj.lig            # must match the NAME in the params file exactly
    motifs, nproblematic, ntries = [], 0, 0
    while len(motifs) < int(m["n_polar_motifs"]) and ntries < int(m["max_attempts"]):
        problem_res, ngly, cur, used = 0, 0, [], []
        for s in spec:
            donor = bool(s.get("donor", False))
            hyb = HYB.get(str(s.get("hybridization", "sp2")), "y")
            if donor:
                res = random.choice(ACCEPTORS)
            else:
                res = random.choice(Q_POLAR if s.get("charged_residues", True) else POLAR)
            used.append(res)
            ngly += res == "GLY"
            if res != "GLY" and res in PROBLEMATIC:
                problem_res += 1
            cur.append(_hbond_block(res, lig, list(s["atoms"]), donor, hyb))
        uniq_ok = len(set(used)) >= len(spec) - 3
        if problem_res == 1 and nproblematic <= int(m["allowed_problematic_motifs"]) and ngly <= 2 and uniq_ok and cur not in motifs:
            nproblematic += 1
            motifs.append(cur)
        elif problem_res == 0 and ngly <= 2 and uniq_ok and cur not in motifs:
            motifs.append(cur)
        ntries += 1
    return motifs


def _write_cst(path: Path, blocks: list[list[str]]) -> None:
    with open(path, "w") as of:
        for b in blocks:
            of.write("CST::BEGIN\n" + "\n".join(b) + "\nCST::END\n")


# ------------------------------------------------------------------ geometry (from workflow.py)
def _angle(p1, p2, p3):
    a, b, c = map(np.array, (p1, p2, p3))
    ba, bc = a - b, c - b
    return float(np.degrees(np.arccos(np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc)))))


def _dihedral(i1, i2, i3, i4):
    p0, p1, p2, p3 = map(np.array, (i1, i2, i3, i4))
    b0, b1, b2 = -1.0 * (p1 - p0), p2 - p1, p3 - p2
    b0xb1, b1xb2 = np.cross(b0, b1), np.cross(b2, b1)
    y = np.dot(np.cross(b0xb1, b1xb2), b1) * (1.0 / np.linalg.norm(b1))
    return float(np.degrees(np.arctan2(y, np.dot(b0xb1, b1xb2))))


def _cluster_paths(proj: Project) -> list[tuple[str, int]]:
    m = proj.config["motifs"]
    cdir = proj.stage("cluster", False)
    names = {f["name"]: f"Fragment_{i}" for i, f in enumerate(proj.fragments(), 1)}
    hb = {names.get(n, n) for n in proj.config.get("hbond_fragments") or []}
    np_frags = [names.get(n, n) for n in (m.get("nonpolar_fragments") or [])] or \
               [p.name for p in cdir.iterdir() if p.is_dir() and p.name.startswith("Fragment_") and p.name not in hb]
    paths = []
    for frag in np_frags:
        for c in (cdir / frag).glob("*_clustered/*.pdb"):
            if c.name.split("_")[1] in m["nonpolar_residues"]:
                paths.append(str(c))
    for frag in hb:   # H-bonding clusters are also offered as explicit contacts, as in workflow.py
        paths += [str(c) for c in (cdir / frag / "hbond_clusts").glob("*.pdb")]
    out = [(p, int(Path(p).stem.split("_")[3])) for p in paths]
    return sorted(out, key=lambda x: -x[1])


def nonpolar_contacts(proj: Project, outdir: Path):
    from pyrosetta import AtomID, Pose, ScoreFunction, generate_nonstandard_residue_set, init, pose_from_file
    from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep
    init("-pack_missing_sidechains False -mute all")
    m = proj.config["motifs"]
    sf = ScoreFunction()
    sf.set_weight(fa_atr, 1)
    sf.set_weight(fa_rep, .55)
    params = [str(proj.params())]
    outdir.mkdir(parents=True, exist_ok=True)
    weighted, named, pdbout = [], [], []
    for clust, freq in _cluster_paths(proj)[:int(m["n_top_clusters"])]:
        p = Pose()
        generate_nonstandard_residue_set(p, params)
        pose_from_file(p, clust)
        lig_i = p.total_residue()
        scores = []
        for i in range(1, lig_i):
            cp = Pose()
            cp.append_residue_by_jump(p.residue(i).clone(), 1)
            cp.append_residue_by_jump(p.residue(lig_i).clone(), 1)
            if cp.residue(1).name()[:3] == cp.residue(2).name()[:3]:
                continue
            scores.append((i, sf(cp)))
        if not scores:
            continue
        best_i, best = sorted(scores, key=lambda x: x[1])[0]
        if best > float(m["single_contact_threshold"]):
            continue
        cp = Pose()
        cp.append_residue_by_jump(p.residue(best_i).clone(), 1)
        cp.append_residue_by_jump(p.residue(lig_i).clone(), 1)
        pdbpath = outdir / f"npmotif_{len(pdbout) + 1}.pdb"
        res_name, lig_name = cp.residue(1).name()[:3], cp.residue(2).name()[:3]
        la = {cp.residue(2).atom_name(k).strip(): (*cp.residue(2).atom(k).xyz(), k) for k in range(1, cp.residue(2).natoms())}
        ra = {cp.residue(1).atom_name(j).strip(): (*cp.residue(1).atom(j).xyz(), j) for j in range(1, cp.residue(1).natoms())}
        d = sorted(((k1, k2, math.dist(v1[:3], v2[:3]), v1[3], v2[3]) for k1, v1 in la.items() for k2, v2 in ra.items()
                    if "H" not in k1 and "H" not in k2), key=lambda x: x[2])
        pick = next((x for x in d[:5] if x[1] not in ("N", "C", "CA", "O")), d[min(4, len(d) - 1)])
        lig_a1, res_a1, dAB, il, ir = pick

        def chain(resi, start):
            nb1 = [a for a in range(1, cp.residue(resi).natoms()) if a != start and
                   cp.conformation().is_bonded(AtomID(start, resi), AtomID(a, resi))]
            for x in nb1:
                for a2 in range(1, cp.residue(resi).natoms()):
                    if a2 not in (x, start) and cp.conformation().is_bonded(AtomID(x, resi), AtomID(a2, resi)) \
                            and "H" not in cp.residue(resi).atom_name(a2).strip():
                        return x, a2
            return None
        lc, rc = chain(2, il), chain(1, ir)
        if lc is None or rc is None:
            continue
        lig_atoms = [lig_a1, cp.residue(2).atom_name(lc[0]).strip(), cp.residue(2).atom_name(lc[1]).strip()]
        res_atoms = ["O" if a == "OXT" else a for a in
                     [res_a1, cp.residue(1).atom_name(rc[0]).strip(), cp.residue(1).atom_name(rc[1]).strip()]]
        if len(set(res_atoms)) < 3 or len(set(lig_atoms)) < 3:
            continue
        cp.dump_pdb(str(pdbpath))
        pdbout.append(str(pdbpath))
        R = [tuple(cp.residue(1).atom(a).xyz()) for a in (ir, *rc)]
        L = [tuple(cp.residue(2).atom(a).xyz()) for a in (il, *lc)]
        block = nonpolar_block(dAB, _angle(L[1], L[0], R[0]), _angle(L[0], R[0], R[1]), _dihedral(L[2], L[1], L[0], R[0]),
                               _dihedral(L[1], L[0], R[0], R[1]), _dihedral(L[0], R[0], R[1], R[2]),
                               res_name, lig_name, lig_atoms, res_atoms,
                               **_apply("nonpolar", dict(distance_tolerance=0.25, angle_A_tolerance=15, angle_B_tolerance=15,
                                                         torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                                         torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                         distance_constraint_sample_number=2)))
        weighted += [block] * freq            # population-weighted sampling, as in workflow.py
        named.append((str(pdbpath), block))
    # pairs of compatible nonpolar contacts
    pairs = []
    if len(pdbout) > 1:
        poses = {}
        for f in pdbout:
            q = Pose()
            generate_nonstandard_residue_set(q, params)
            pose_from_file(q, f)
            poses[f] = q
        sc = []
        for i, a in enumerate(pdbout):
            for b in pdbout[i + 1:]:
                cp = Pose()
                cp.append_residue_by_jump(poses[a].residue(1).clone(), 1)
                cp.append_residue_by_jump(poses[b].residue(1).clone(), 1)
                cp.append_residue_by_jump(poses[a].residue(2).clone(), 1)
                sc.append((a, b, sf(cp)))
        blk = dict(named)
        pairs = [(blk[a], blk[b]) for a, b, s in sorted(sc, key=lambda x: x[2]) if s <= float(m["double_contact_threshold"])]
    return weighted, pairs, len(pdbout)


def run(proj: Project, args) -> None:
    m = proj.config["motifs"]
    random.seed(int(m.get("seed", 1)))
    out = proj.stage("motifs")
    fast = bool(getattr(args, "fast", False))
    cons = (m.get("fast_preset") if fast else m.get("constraints")) or {}
    set_overrides(cons.get("hbond"), cons.get("nonpolar"))
    mode = "fast preset" if fast else ("custom overrides" if any(cons.get(k) for k in ("hbond", "nonpolar")) else "workflow.py defaults")
    print(f"[motifs] constraint geometry: {mode}")
    (out / "constraint_settings.json").write_text(__import__("json").dumps({"mode": mode, "overrides": cons}, indent=1))
    polar = polar_motifs(proj)
    pdir = out / "polar"
    pdir.mkdir(exist_ok=True)
    for old in pdir.glob("*.cst"):
        old.unlink()
    for i, mo in enumerate(polar, 1):
        _write_cst(pdir / f"cst_{i}.cst", mo)
    print(f"[motifs] {len(polar)} polar motifs in {pdir}")
    weighted, pairs, n_np = nonpolar_contacts(proj, out / "npmotif_pdbs")
    print(f"[motifs] {n_np} nonpolar contacts accepted, {len(pairs)} compatible pairs")
    hybrids = {int(k): int(v) for k, v in m["hybrids"].items()}
    bad = [k for k in hybrids if k not in (1, 2)]
    if bad:
        raise ClaireError(f"motifs.hybrids keys must be 1 or 2 (number of nonpolar contacts per motif); got {bad}.")
    sets = {}
    for n, per in hybrids.items():
        pool = weighted if n == 1 else pairs
        d = out / f"{n}np"
        d.mkdir(exist_ok=True)
        for old in d.glob("*.cst"):
            old.unlink()
        if not pool:
            print(f"[motifs] no {'contacts' if n == 1 else 'contact pairs'} available for {n}np motifs")
            continue
        c = 0
        for mo in polar:
            chosen = []
            for _ in range(per * 3):
                if len(chosen) == per:
                    break
                r = random.choice(pool)
                if r not in chosen:
                    chosen.append(r)
            for r in chosen:
                c += 1
                extra = [r] if n == 1 else list(r)
                _write_cst(d / f"hybrid_{c}.cst", extra + mo)
        sets[f"{n}np"] = c
        print(f"[motifs] {c} motifs with {n} nonpolar contact(s) in {d}")
    proj.write_manifest("motifs", {"polar": len(polar), "sets": sets, "constraints": mode, "n_outputs": sum(sets.values())})
    print("Next: `claire match` (scaffolds from scaffolds.dir).")
