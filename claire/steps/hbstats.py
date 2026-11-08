"""Statistical hydrogen-bond geometry for motif constraints (motifs.hbond_geometry: statistical).

For each ligand polar atom in `polar_atoms`, contacts mined for the fragments containing it
(03_contacts/residue_contact_fuzzballs) are scanned for heavy-atom hydrogen bonds (2.5-3.3 A, complementary
donor/acceptor roles). Per (residue type, residue polar atom), the six matcher degrees of freedom are measured between
heavy atoms (the contact structures carry no protein hydrogens):

    distanceAB  L1...R1            angle_A  L2-L1...R1        angle_B  L1...R1-R2
    torsion_A   L3-L2-L1...R1      torsion_AB  L2-L1...R1-R2  torsion_B  L1...R1-R2-R3

Ideal values are medians (angles, distance) or circular means (torsions); tolerances come from the spread
(1.5 x robust SD, clipped to sensible ranges); torsions with no clear preference stay free. Residue types are then
sampled in proportion to their observed counts.
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from ..project import ClaireError, Project

# residue polar atoms: (atom, bonded, next, role); symmetric atoms are folded onto one name
RES_POLAR = {
    "SER": [("OG", "CB", "CA", "both")], "THR": [("OG1", "CB", "CA", "both")], "TYR": [("OH", "CZ", "CE1", "both")],
    "TRP": [("NE1", "CD1", "CG", "don")], "LYS": [("NZ", "CE", "CD", "don")],
    "ARG": [("NE", "CD", "CG", "don"), ("NH1", "CZ", "NE", "don"), ("NH2", "CZ", "NE", "don")],
    "ASN": [("OD1", "CG", "CB", "acc"), ("ND2", "CG", "CB", "don")],
    "GLN": [("OE1", "CD", "CG", "acc"), ("NE2", "CD", "CG", "don")],
    "HIS": [("ND1", "CG", "CB", "both"), ("NE2", "CD2", "CG", "both")],
    "ASP": [("OD1", "CG", "CB", "acc"), ("OD2", "CG", "CB", "acc")],
    "GLU": [("OE1", "CD", "CG", "acc"), ("OE2", "CD", "CG", "acc")],
    "GLY": [("N", "CA", "C", "don"), ("O", "C", "CA", "acc")],
}
FOLD = {("ARG", "NH2"): "NH1", ("ASP", "OD2"): "OD1", ("GLU", "OE2"): "OE1"}
CHARGED = {"LYS", "ARG", "ASP", "GLU"}


def _ang(a, b, c):
    v1, v2 = a - b, c - b
    return math.degrees(math.acos(max(-1.0, min(1.0, float(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)))))))


def _dih(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return math.degrees(math.atan2(float(np.cross(b1, v) @ w), float(v @ w)))


def _read(pdb: Path, lig: str):
    ligand, residues = {}, defaultdict(dict)
    for l in pdb.read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        name, rn = l[12:16].strip(), l[17:20].strip()
        if rn.upper() == lig.upper():
            ligand[name] = xyz
        else:
            residues[(rn, l[21], l[22:27])][name] = xyz
    return ligand, residues


def ligand_triplet(proj: Project, spec: dict) -> tuple[list[str], bool]:
    """(heavy-atom triplet L1, L2, L3, ligand_is_donor) for a polar_atoms entry."""
    from .autofrag import chemistry
    atoms = list(spec["atoms"])
    donor = bool(spec.get("donor", False))
    if not donor:
        return atoms[:3], False
    _, heavy, el, nb, _, _ = chemistry(proj)
    l1, l2 = atoms[1], atoms[2]
    l3 = next((x for x in nb[l2] if x != l1 and el.get(x) != "H"), None)
    if l3 is None:
        raise ClaireError(f"Cannot find a third heavy atom for donor {atoms}.")
    return [l1, l2, l3], True


def _fit(vals: np.ndarray, kind: str):
    """(ideal, tolerance, samples) for one degree of freedom."""
    if kind == "torsion":
        r = np.radians(vals)
        c, s = np.cos(r).mean(), np.sin(r).mean()
        R = math.hypot(c, s)
        if R < 0.3:
            return 0.0, 360.0, 9
        sd = math.degrees(math.sqrt(-2 * math.log(R)))
        tol = float(np.clip(1.5 * sd, 15, 60))
        return math.degrees(math.atan2(s, c)), tol, int(min(4, math.ceil(tol / 15)))
    med = float(np.median(vals))
    mad = float(np.median(np.abs(vals - med))) * 1.4826
    if kind == "distance":
        tol = float(np.clip(1.5 * mad, 0.15, 0.5))
        return med, tol, int(min(3, math.ceil(tol / 0.15)))
    tol = float(np.clip(1.5 * mad, 10, 35))
    return med, tol, int(min(3, math.ceil(tol / 10)))


def analyse(proj: Project, spec: dict, min_n: int) -> tuple[list[dict], list[dict]]:
    """Fitted options for one polar atom: [{res, ratom, n, params...}] and all raw observations."""
    (L1, L2, L3), lig_donor = ligand_triplet(proj, spec)
    names = {f["name"]: (f"Fragment_{i}", f["atoms"]) for i, f in enumerate(proj.fragments(), 1)}
    frags = [fid for fid, atoms in names.values() if {L1, L2}.issubset(atoms)]
    base = proj.stage("contacts", False) / "residue_contact_fuzzballs"
    files = []
    for fid in frags:
        files += sorted((base / f"{proj.lig}_{fid}_residue_contacts").glob("*.pdb"))
    obs = defaultdict(list)
    raw = []
    for f in files:
        ligand, residues = _read(f, proj.lig)
        if not all(a in ligand for a in (L1, L2, L3)):
            continue
        Lp = [ligand[a] for a in (L1, L2, L3)]
        for (rn, ch, num), at in residues.items():
            for ra, rb, rc, role in RES_POLAR.get(rn, []):
                if not all(x in at for x in (ra, rb, rc)):
                    continue
                if (lig_donor and role == "don") or (not lig_donor and role == "acc"):
                    continue
                d = float(np.linalg.norm(Lp[0] - at[ra]))
                if not 2.5 <= d <= 3.3:
                    continue
                R = [at[ra], at[rb], at[rc]]
                dof = {"distance": d, "angle_A": _ang(Lp[1], Lp[0], R[0]), "angle_B": _ang(Lp[0], R[0], R[1]),
                       "torsion_A": _dih(Lp[2], Lp[1], Lp[0], R[0]), "torsion_AB": _dih(Lp[1], Lp[0], R[0], R[1]),
                       "torsion_B": _dih(Lp[0], R[0], R[1], R[2])}
                key = (rn, FOLD.get((rn, ra), ra))
                obs[key].append(dof)
                raw.append({"ligand_atom": L1, "residue": rn, "residue_atom": key[1], "source": f.name, **dof})
    options = []
    for (rn, ra), rows in sorted(obs.items(), key=lambda kv: -len(kv[1])):
        if len(rows) < min_n:
            continue
        o = {"res": rn, "ratom": ra, "rtrip": next([a, b, c] for a, b, c, _ in RES_POLAR[rn] if FOLD.get((rn, a), a) == ra),
             "ltrip": [L1, L2, L3], "n": len(rows)}
        for k, kind in (("distance", "distance"), ("angle_A", "angle"), ("angle_B", "angle"),
                        ("torsion_A", "torsion"), ("torsion_AB", "torsion"), ("torsion_B", "torsion")):
            o[k] = _fit(np.array([r[k] for r in rows]), kind)
        options.append(o)
    return options, raw


def block(o: dict, lig: str, apply) -> list[str]:
    """Matcher constraint block with explicit heavy-atom names on both sides (periodicity 360)."""
    kw = {"distance_tolerance": o["distance"][1], "angle_A_tolerance": o["angle_A"][1],
          "angle_B_tolerance": o["angle_B"][1], "torsion_A_tolerance": o["torsion_A"][1],
          "torsion_AB_tolerance": o["torsion_AB"][1], "torsion_B_tolerance": o["torsion_B"][1],
          "distance_constraint_sample_number": o["distance"][2], "angle_constraint_sample_number": max(o["angle_A"][2], o["angle_B"][2]),
          "torsion_A_samples": o["torsion_A"][2], "torsion_AB_samples": o["torsion_AB"][2], "torsion_B_samples": o["torsion_B"][2],
          "torsion_constraint_sample_number": o["torsion_A"][2]}
    kw = apply("hbond", kw)
    line = "  CONSTRAINT:: {0:>10}: {1:7.2f} {2:6.2f} {3:6.2f}  {4:>6}  {5:3}"
    return [
        "  TEMPLATE::   ATOM_MAP: 1 atom_name: {}".format(" ".join(o["ltrip"])),
        "  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n".format(lig),
        "  TEMPLATE::   ATOM_MAP: 2 atom_name: {}".format(" ".join(o["rtrip"])),
        "  TEMPLATE::   ATOM_MAP: 2 residue3: {}\n".format(o["res"]),
        line.format("distanceAB", o["distance"][0], kw["distance_tolerance"], 100, "0", kw["distance_constraint_sample_number"]),
        line.format("angle_A", o["angle_A"][0], kw["angle_A_tolerance"], 100, "360.00", kw["angle_constraint_sample_number"]),
        line.format("angle_B", o["angle_B"][0], kw["angle_B_tolerance"], 100, "360.00", kw["angle_constraint_sample_number"]),
        line.format("torsion_A", o["torsion_A"][0], kw["torsion_A_tolerance"], 100, "360.00", kw["torsion_A_samples"]),
        line.format("torsion_B", o["torsion_B"][0], kw["torsion_B_tolerance"], 100, "360.00", kw["torsion_B_samples"]),
        line.format("torsion_AB", o["torsion_AB"][0], kw["torsion_AB_tolerance"], 100, "360.00", kw["torsion_AB_samples"]),
    ]


def write_report(path: Path, rows: list[dict]) -> None:
    keys = ["ligand_atom", "residue", "residue_atom", "n", "weight", "distance", "angle_A", "angle_B", "torsion_A",
            "torsion_AB", "torsion_B", "status"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
