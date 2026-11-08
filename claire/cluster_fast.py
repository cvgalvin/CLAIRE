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


def main(argv=None):
    argv = argv or sys.argv[1:]
    rfb, ligand_path = argv[0], argv[1]
    lines, residues, rows = read_fuzzball(rfb)
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
