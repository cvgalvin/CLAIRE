"""Fast drop-in replacement for the original spatial_cluster.py (same arguments, clusters, files and file contents).

usage: python -m claire.cluster_fast <residue_contact_fuzzball.pdb> <ligand.pdb>

Algorithm (unchanged): greedy leader clustering in the ligand frame. Residues are visited in file order; each residue not
yet clustered becomes a centre, and every later unclustered residue whose RMSD to the centre over shared side-chain
heavy atoms (no superposition) is <= 0.5 A joins it. Only clusters with at least two members are written.

Speed-ups (results identical):
  * O(1) bookkeeping of clustered residues (the original list lookup made the scan O(n^3));
  * vectorized RMSD of a centre against all candidates;
  * exact pre-filter: RMSD <= 0.5 A over identical atom sets implies centroid distance <= 0.5 A, so candidates for
    residues with the full atom set come from a KD-tree on centroids; residues with missing atoms are compared
    against all candidates over their shared atoms, as in the original;
  * residues with duplicate atom names (alternate locations) use the original pairwise routine.
"""
from __future__ import annotations

import os
import sys
from collections import OrderedDict

import numpy as np

BACKBONE = {"OXT", "C", "CA", "O", "N"}
CUTOFF = 0.5


def _orig_line(line: str) -> str:
    """Line exactly as the original df_to_pdb wrote it (fields re-joined, charge columns dropped)."""
    cuts = [(0, 6), (6, 11), (11, 16), (16, 17), (17, 20), (20, 22), (22, 29), (29, 31), (31, 39), (39, 47), (47, 54),
            (54, 60), (60, 66), (66, 76), (76, 78)]
    return "".join(line[a:b] for a, b in cuts) + "\n"


def _orig_lig_line(line: str) -> str:
    cuts = [(6, 11), (11, 16), (16, 17), (17, 20)]
    rest = [(22, 29), (29, 31), (31, 39), (39, 47), (47, 54), (54, 60), (60, 66), (66, 76), (76, 78)]
    return "ATOM  " + "".join(line[a:b] for a, b in cuts) + " X" + "".join(line[a:b] for a, b in rest) + "\n"


def read_fuzzball(path):
    lines = [l for l in open(path) if l[0:4] == "ATOM"]
    residues = OrderedDict()          # resnum field -> list of (atom name, x, y, z), deduplicated as in the original
    rows = OrderedDict()              # resnum field -> list of line indices
    for i, line in enumerate(lines):
        key = line[22:29]
        rows.setdefault(key, []).append(i)
        name = line[11:16].split()[0]
        if name[0] != "H" and name not in BACKBONE:
            info = (name, float(line[31:39]), float(line[39:47]), float(line[47:54]))
            lst = residues.setdefault(key, [])
            if info not in lst:
                lst.append(info)
    return lines, residues, rows


def _pair_rmsd_original(r1, r2):
    vs1, vs2 = [], []
    for a, b, c, d in r1:
        for e, f, g, h in r2:
            if a == e:
                vs1.append((b, c, d))
                vs2.append((f, g, h))
    if not vs1:
        return float("nan")
    diff = np.array(vs1) - np.array(vs2)
    return float(np.sqrt((diff * diff).sum() / len(vs1)))


def cluster(residues: "OrderedDict[str, list]") -> "OrderedDict[str, list]":
    keys = list(residues)
    n = len(keys)
    names = sorted({a[0] for r in residues.values() for a in r})
    col = {nm: j for j, nm in enumerate(names)}
    X = np.full((n, len(names), 3), np.nan)
    dup = np.zeros(n, bool)
    for i, k in enumerate(keys):
        seen = set()
        for nm, x, y, z in residues[k]:
            if nm in seen:
                dup[i] = True
                continue
            seen.add(nm)
            X[i, col[nm]] = (x, y, z)
    present = ~np.isnan(X[:, :, 0])
    complete = present.all(axis=1) & ~dup
    tree = None
    if complete.any():
        from scipy.spatial import cKDTree
        cidx = np.where(complete)[0]
        cent = X[cidx].mean(axis=1)
        tree = cKDTree(cent)
    incomplete_idx = np.where(~complete)[0]
    clustered = np.zeros(n, bool)
    out = OrderedDict()
    for i in range(n):
        if clustered[i]:
            continue
        clustered[i] = True
        if dup[i]:
            cand = np.where(~clustered)[0]
            cand = cand[cand > i]
            members = [j for j in cand if _pair_rmsd_original(residues[keys[i]], residues[keys[j]]) <= CUTOFF]
        else:
            if complete[i] and tree is not None:
                near = cidx[tree.query_ball_point(X[i].mean(axis=0), CUTOFF + 1e-9)]
                cand = np.concatenate([near, incomplete_idx])
            else:
                cand = np.arange(n)
            cand = np.unique(cand[(cand > i)])
            cand = cand[~clustered[cand]]
            members = []
            if len(cand):
                regular = cand[~dup[cand]]
                if len(regular):
                    d2 = ((X[regular] - X[i][None]) ** 2).sum(axis=2)               # (m, atoms), NaN if missing
                    shared = ~np.isnan(d2)
                    cnt = shared.sum(axis=1)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        rmsd = np.sqrt(np.where(shared, d2, 0).sum(axis=1) / cnt)
                    members += list(regular[(cnt > 0) & (rmsd <= CUTOFF)])
                for j in cand[dup[cand]]:
                    if _pair_rmsd_original(residues[keys[i]], residues[keys[j]]) <= CUTOFF:
                        members.append(j)
            members = sorted(int(j) for j in members)
        if members:
            clustered[members] = True
            out[keys[i]] = [keys[j] for j in members]
    return out


# ---------------------------------------------------------------------------------------------------------------
# Group clustering (v0.4): each residue is compared on its interacting functional group (plus one anchor atom), in
# every symmetry-equivalent orientation of the fragment and every equivalent naming of the residue's atoms.
# ---------------------------------------------------------------------------------------------------------------
GROUP_CUTOFF = 0.5
GROUPS = {   # residue: (group atoms, anchor first; symmetric renamings, each a set of name pairs exchanged together)
    "ARG": (["CD", "NE", "CZ", "NH1", "NH2"], [(("NH1", "NH2"),)]),
    "LYS": (["CD", "CE", "NZ"], []),
    "ASP": (["CB", "CG", "OD1", "OD2"], [(("OD1", "OD2"),)]),
    "GLU": (["CG", "CD", "OE1", "OE2"], [(("OE1", "OE2"),)]),
    "ASN": (["CB", "CG", "OD1", "ND2"], []),
    "GLN": (["CG", "CD", "OE1", "NE2"], []),
    "SER": (["CA", "CB", "OG"], []),
    "THR": (["CA", "CB", "OG1", "CG2"], []),
    "CYS": (["CA", "CB", "SG"], []),
    "MET": (["CG", "SD", "CE"], []),
    "HIS": (["CB", "CG", "ND1", "CD2", "CE1", "NE2"], []),
    "PHE": (["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"], [(("CD1", "CD2"), ("CE1", "CE2"))]),
    "TYR": (["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"], [(("CD1", "CD2"), ("CE1", "CE2"))]),
    "TRP": (["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"], []),
    "LEU": (["CB", "CG", "CD1", "CD2"], [(("CD1", "CD2"),)]),
    "ILE": (["CB", "CG1", "CG2", "CD1"], []),
    "VAL": (["CA", "CB", "CG1", "CG2"], [(("CG1", "CG2"),)]),
    "ALA": (["N", "CA", "CB"], []),
    "PRO": (["N", "CB", "CG", "CD"], []),
}


def name_variants(atoms, swaps):
    """The residue's symmetric renamings (e.g. Phe ring flip: CD1<->CD2 together with CE1<->CE2)."""
    vs = [list(atoms)]
    for pairs in swaps:
        m = {}
        for a, b in pairs:
            m[a], m[b] = b, a
        vs = vs + [[m.get(x, x) for x in v] for v in vs]
    out = []
    for v in vs:
        if v not in out:
            out.append(v)
    return out


def _kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Qc - Pc @ R.T


def _cc_bond_types(ligand_path, names):
    """Bond types between carbon atoms, from the ligand's params file next to the model (bond orders with RDKit
    aromaticity perception, as for fragment search patterns); None if unavailable."""
    import glob, os
    from rdkit import Chem
    params = sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(ligand_path)), "*.params")))
    if not params:
        return None
    atoms, el, bonds = [], {}, []
    for l in open(params[0]):
        p = l.split()
        if l.startswith("ATOM "):
            atoms.append(p[1])
        elif l.startswith("BOND_TYPE"):
            bonds.append((p[1], p[2], p[3]))
        elif l.startswith("BOND ") and len(p) >= 3:
            bonds.append((p[1], p[2], "1"))
    for l in open(ligand_path):
        if l.startswith(("HETATM", "ATOM")):
            n = l[12:16].strip()
            el[n] = (l[76:78].strip() or n.lstrip("0123456789")[0]).capitalize()
    try:
        m = Chem.RWMol(); idx = {}
        for a in atoms:
            at = Chem.Atom(el[a]); at.SetNoImplicit(True); idx[a] = m.AddAtom(at)
        bt = {"1": Chem.BondType.SINGLE, "2": Chem.BondType.DOUBLE, "3": Chem.BondType.TRIPLE, "4": Chem.BondType.AROMATIC}
        for a, b, o in bonds:
            if m.GetBondBetweenAtoms(idx[a], idx[b]) is None:
                m.AddBond(idx[a], idx[b], bt.get(o, Chem.BondType.SINGLE))
        mol = m.GetMol()
        for a in mol.GetAtoms():                      # formal charges implied by valence, as in autofrag.parent_mol
            if a.GetSymbol() == "O" and a.GetDegree() == 1 and mol.GetBondBetweenAtoms(a.GetIdx(), a.GetNeighbors()[0].GetIdx()).GetBondType() == Chem.BondType.SINGLE \
                    and a.GetNeighbors()[0].GetSymbol() != "H":
                a.SetFormalCharge(-1)
            if a.GetSymbol() == "N" and not any(b.GetBondType() == Chem.BondType.AROMATIC for b in a.GetBonds()) \
                    and sum(b.GetBondTypeAsDouble() for b in a.GetBonds()) == 4:
                a.SetFormalCharge(1)
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    out = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if a in idx and b in idx and el.get(a) == "C" and el.get(b) == "C":
                bd = mol.GetBondBetweenAtoms(idx[a], idx[b])
                if bd is not None:
                    out[frozenset((a, b))] = bd.GetBondType()
    return out


def fragment_symmetry(ligand_path, fragment_atoms):
    """Rigid transforms (R, t) mapping the fragment onto itself, one per automorphism of its heavy-atom graph
    (12 for benzene, 6 for guanidine, 2 for a carboxylate or a phenol ring flip, 1 for an asymmetric fragment).
    Carbon-carbon bonds are typed (from the params' bond orders with aromaticity perception, else from bond length),
    so a C=C is not equivalent to a C-C (cyclohexenone has no mirror symmetry); bonds to heteroatoms are left untyped,
    so resonance-equivalent atoms stay interchangeable (carboxylate and sulfonate oxygens, guanidinium nitrogens)."""
    names, xyz, els = [], [], []
    for l in open(ligand_path):
        if l.startswith(("HETATM", "ATOM")):
            n = l[12:16].strip()
            if n in fragment_atoms:
                names.append(n)
                xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                els.append((l[76:78].strip() or n.lstrip("0123456789")[0]).capitalize())
    X = np.array(xyz)
    if len(X) < 3:
        return [(np.eye(3), np.zeros(3))]
    from rdkit import Chem
    typed = _cc_bond_types(ligand_path, names)
    rw = Chem.RWMol()
    for e in els:
        rw.AddAtom(Chem.Atom(e))
    for a in range(len(X)):
        for b in range(a + 1, len(X)):
            d = np.linalg.norm(X[a] - X[b])
            if d < 1.75:
                t = Chem.BondType.SINGLE
                if els[a] == "C" and els[b] == "C":
                    if typed is not None and frozenset((names[a], names[b])) in typed:
                        t = typed[frozenset((names[a], names[b]))]
                    elif typed is None:
                        t = (Chem.BondType.TRIPLE if d < 1.25 else Chem.BondType.DOUBLE if d < 1.37
                             else Chem.BondType.AROMATIC if d < 1.44 else Chem.BondType.SINGLE)
                rw.AddBond(a, b, t)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    ops = []
    for p in m.GetSubstructMatches(m, uniquify=False, useChirality=False, maxMatches=500):
        R, t = _kabsch(X, X[list(p)])
        if np.sqrt((((X @ R.T + t) - X[list(p)]) ** 2).sum(1).mean()) < 0.3:     # proper rigid symmetries only
            ops.append((R, t))
    return ops or [(np.eye(3), np.zeros(3))]


def cluster_groups(lines, rows, resname, ops, cutoff=GROUP_CUTOFF):
    """Greedy clustering (same order and output as cluster()) on the residue's functional group; the RMSD between two
    contacts is the minimum over the fragment's symmetry operations and the residue's symmetric atom namings.
    Contacts lacking any group atom stay unclustered."""
    from scipy.spatial import cKDTree
    keys = list(rows)
    atoms, swaps = GROUPS[resname]
    variants = name_variants(atoms, swaps)
    coords, ok = [], []
    for k in keys:
        A = {}
        for i in rows[k]:
            A.setdefault(lines[i][12:16].strip(), np.array([float(lines[i][30:38]), float(lines[i][38:46]),
                                                            float(lines[i][46:54])]))
        if all(a in A for a in atoms):
            coords.append([[A[a] for a in v] for v in variants])
            ok.append(True)
        else:
            coords.append(None)
            ok.append(False)
    idx = [i for i in range(len(keys)) if ok[i]]
    out = OrderedDict()
    if len(idx) < 2:
        return out
    X = np.array([coords[i] for i in idx])                                        # (n, variants, atoms, 3)
    imgs = np.einsum("nvad,ode->novae", X, np.array([R.T for R, t in ops])) + np.array([t for R, t in ops])[None, :, None, None, :]
    ref = X[:, 0]
    cent = imgs.mean(axis=3).reshape(len(idx), -1, 3)
    owner = np.repeat(np.arange(len(idx)), cent.shape[1])
    tree = cKDTree(cent.reshape(-1, 3))
    done = np.zeros(len(idx), bool)
    for a in range(len(idx)):
        if done[a]:
            continue
        done[a] = True
        cand = np.unique(owner[tree.query_ball_point(ref[a].mean(0), cutoff + 1e-9)])
        cand = cand[~done[cand]]
        members = [int(b) for b in cand
                   if np.sqrt(((imgs[b] - ref[a]) ** 2).sum(-1).mean(-1)).min() <= cutoff]
        if members:
            done[members] = True
            out[keys[idx[a]]] = [keys[idx[b]] for b in sorted(members)]
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="cluster one per-residue contact file")
    ap.add_argument("fuzzball"); ap.add_argument("ligand")
    ap.add_argument("--method", choices=["group", "sidechain"], default="sidechain",
                    help="group: functional group + symmetry (v0.4 default via claire cluster); sidechain: original")
    ap.add_argument("--cutoff", type=float, default=GROUP_CUTOFF)
    ap.add_argument("--fragment-atoms", default="", help="comma-separated fragment heavy atoms (group method)")
    a = ap.parse_args(argv if argv is not None else sys.argv[1:])
    rfb, ligand_path = a.fuzzball, a.ligand
    lines, residues, rows = read_fuzzball(rfb)
    resname = lines[0][17:20] if lines else ""
    if a.method == "group" and resname in GROUPS:
        ops = fragment_symmetry(ligand_path, set(x for x in a.fragment_atoms.split(",") if x))
        clusters = cluster_groups(lines, rows, resname, ops, a.cutoff)
    else:
        clusters = cluster(residues)
    pops = sorted(((k, len(v) + 1) for k, v in clusters.items()), key=lambda x: -x[1])
    try:                                   # population plot, as in the original
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages(rfb.split(".")[0] + "_clustered.pdf") as pdf:
            plt.bar([i + 1 for i in range(len(pops))], [p for _, p in pops])
            plt.xlabel("Cluster ID")
            plt.ylabel("Number of Members")
            plt.title("Cluster Population Distribution")
            pdf.savefig()
            plt.clf()
    except Exception:
        pass
    if not pops:
        print("there are 0 clusters")
        return
    outdir = rfb.split(".")[0] + "_clustered"
    os.makedirs(outdir, exist_ok=True)
    lig = []
    for l in open(ligand_path):
        if l[0:6] == "HETATM":
            lig.append(_orig_lig_line(l))
        else:
            break
    for key, pop in pops:
        text = []
        for rk in clusters[key] + [key]:
            idx = rows[rk]
            text += [_orig_line(lines[t]) for t in range(min(idx), max(idx) + 1)]
        resname = lines[rows[key][0]][17:20]
        with open(os.path.join(outdir, f"{key.strip()}_{resname}_cluster_{pop}.pdb"), "w") as fh:
            fh.write("".join(text + lig))
    print(f"there are {sum(p for _, p in pops)} contacts after clustering ({len(pops)} clusters)")


if __name__ == "__main__":
    main()
