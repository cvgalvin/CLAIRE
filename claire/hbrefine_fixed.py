"""HBRefine, fixed version (the published algorithm is in claire/legacy/hbrefine.py).

usage: python -m claire.hbrefine_fixed LIG.params design.pdb output_dir [max_energy_increase_REU] [diagnose] [no-minimize] [no-backbone]

Same aim and output names as the original (<name>_unrefined.pdb, <name>_<n>mutations_fixed/_unsats/_bothproblems.pdb),
with these changes:
  * binding-site residues: residues within 10 A of the ligand with SASA < 20 A^2, filtered in one correct pass
  * ligand polar atoms: every ligand heavy atom Rosetta annotates as an H-bond acceptor or donor (any atom type)
  * burial: only *buried* ligand polar atoms count as unsatisfied (atom SASA in the complex < 25% of that in the
    free ligand); solvent-exposed ones are left to water
  * H-bond mutations: candidate residue types depend on what the unsatisfied atoms need (ligand acceptors -> donor
    side chains, ligand donors -> acceptor side chains), and a mutation is accepted only if it reduces the number of
    buried unsatisfied ligand polar atoms
  * "useless" buried polar residues: a polar side chain (incl. His, Tyr, Trp) in the binding site whose buried polar
    atoms form no H-bond with anything (ligand, protein); residues H-bonding the protein are left alone
  * every accepted mutation may raise the total energy by at most max_energy_increase REU (default 1.0)
  * score function: beta_nov16 (as in enzyme design and FastDesign)
"""
from __future__ import annotations

import math
import os
import shutil
import sys

import numpy as np

LIG_BURIED_REL = 0.25        # ligand polar atom is buried if SASA(complex) < 25% of SASA(free ligand)
PROT_BURIED_ABS = 2.0        # protein polar atom is buried if SASA(complex) < 2 A^2
SHELL_CUTOFF = 10.0          # A, binding-site neighbourhood
SHELL_SASA = 20.0            # A^2, residue SASA cutoff for binding-site residues
DONOR_RES = ["S", "T", "Y", "W", "N", "Q", "H"]          # partners for ligand acceptors
ACCEPTOR_RES = ["D", "E", "N", "Q", "H", "S", "T", "Y"]  # partners for ligand donors
POLAR_3 = {"SER", "THR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS", "TYR", "TRP"}
NONPOLAR = ["L", "I", "V", "F", "M"]
# size-matched nonpolar analogs for the control mutation (same bulk, no hydrogen bond)
NP_ANALOG = {"S": "A", "T": "V", "N": "L", "D": "L", "Q": "M", "E": "M", "K": "M", "R": "M", "H": "F", "Y": "F", "W": "F"}
RADII = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}


def _sphere(n=100):
    i = np.arange(n) + 0.5
    phi, th = np.arccos(1 - 2 * i / n), np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


SPHERE = _sphere()


def _sasa(xyz: np.ndarray, rad: np.ndarray, idx, probe=1.4) -> np.ndarray:
    out = []
    R = rad + probe
    for i in idx:
        pts = xyz[i] + SPHERE * R[i]
        d = np.linalg.norm(pts[:, None] - xyz[None], axis=2)
        d[:, i] = 1e9
        out.append(4 * math.pi * R[i] ** 2 * np.mean((d >= R[None]).all(1)))
    return np.array(out)


def _heavy(p, resids):
    """(resid, atom_index, xyz, radius) for heavy atoms of the given residues."""
    rows = []
    for r in resids:
        res = p.residue(r)
        for a in range(1, res.nheavyatoms() + 1):
            el = res.atom_type(a).element().upper()
            v = res.xyz(a)
            rows.append((r, a, (v.x, v.y, v.z), RADII.get(el, 1.8)))
    return rows


def analyse(p, lig: int, shell: list[int]):
    """H-bonds, buried unsatisfied ligand polar atoms, and protein atoms that make any H-bond."""
    from pyrosetta import rosetta
    hbs = rosetta.core.scoring.hbonds.HBondSet()
    p.update_residue_neighbors()
    rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hbs)
    lig_hb_atoms, res_hb_lig, satisfied = set(), set(), set()
    for k in range(1, hbs.nhbonds() + 1):
        hb = hbs.hbond(k)
        dres, ares = int(hb.don_res()), int(hb.acc_res())
        datm = int(p.residue(dres).first_adjacent_heavy_atom(int(hb.don_hatm())))
        aatm = int(hb.acc_atm())
        satisfied.add((dres, datm))
        satisfied.add((ares, aatm))
        if ares == lig and dres != lig:
            lig_hb_atoms.add(aatm)
            res_hb_lig.add(dres)
        elif dres == lig and ares != lig:
            lig_hb_atoms.add(datm)
            res_hb_lig.add(ares)
    L = p.residue(lig)
    polar = [a for a in range(1, L.nheavyatoms() + 1) if L.heavyatom_is_an_acceptor(a) or L.heavyatom_has_polar_hydrogens(a)]
    # burial of ligand polar atoms: complex (ligand + binding-site neighbourhood) vs free ligand
    rows_l = _heavy(p, [lig])
    env = _heavy(p, shell)
    xyz = np.array([r[2] for r in rows_l + env])
    rad = np.array([r[3] for r in rows_l + env])
    pos = {r[1]: i for i, r in enumerate(rows_l)}
    idx = [pos[a] for a in polar]
    s_cx = _sasa(xyz, rad, idx)
    s_free = _sasa(xyz[:len(rows_l)], rad[:len(rows_l)], idx)
    buried = [a for a, c, f in zip(polar, s_cx, s_free) if f <= 0 or c / f < LIG_BURIED_REL]
    unsat = [a for a in buried if a not in lig_hb_atoms]
    roles = {a: ("acc" if L.heavyatom_is_an_acceptor(a) else "") + ("don" if L.heavyatom_has_polar_hydrogens(a) else "")
             for a in unsat}
    return {"unsat": unsat, "roles": roles, "res_hb_lig": res_hb_lig, "satisfied": satisfied}


def useless_polar(p, lig: int, shell: list[int], info) -> list[int]:
    """Binding-site polar side chains whose buried polar atoms form no H-bond at all."""
    rows = _heavy(p, shell + [lig])
    xyz = np.array([r[2] for r in rows])
    rad = np.array([r[3] for r in rows])
    out = []
    for r in shell:
        res = p.residue(r)
        if res.name3() not in POLAR_3 or r in info["res_hb_lig"]:
            continue
        sc_polar = [a for a in range(res.first_sidechain_atom(), res.nheavyatoms() + 1)
                    if res.heavyatom_is_an_acceptor(a) or res.heavyatom_has_polar_hydrogens(a)]
        if not sc_polar:
            continue
        idx = [i for i, row in enumerate(rows) if row[0] == r and row[1] in sc_polar]
        buried = [rows[i][1] for i, s in zip(idx, _sasa(xyz, rad, idx)) if s < PROT_BURIED_ABS]
        if buried and not any((r, a) in info["satisfied"] for a in sc_polar):
            out.append(r)
    return out


def main(argv=None):
    argv = argv or sys.argv[1:]
    params, pdb, outdir = argv[0], argv[1], argv[2]
    cap = float(argv[3]) if len(argv) > 3 else 1.0
    flags = set(argv[4:])
    diagnose = "diagnose" in flags              # record every H-bond-forming candidate, not only those within the cap
    minimize = "no-minimize" not in flags       # minimize the binding site after each H-bond-forming mutation
    min_backbone = "no-backbone" not in flags
    os.makedirs(outdir, exist_ok=True)
    stem = os.path.basename(pdb)[:-4] if pdb.endswith(".pdb") else os.path.basename(pdb)

    import pyrosetta
    from pyrosetta import MoveMap, Pose, generate_nonstandard_residue_set, init, pose_from_file, rosetta
    from pyrosetta.toolbox import mutate_residue
    init("-load_PDB_components False -ex1 -ex2 -beta_nov16 -mute all")
    sf = pyrosetta.create_score_function("beta_nov16")

    p = Pose()
    generate_nonstandard_residue_set(p, [params])
    pose_from_file(p, pdb)
    sf(p)
    ligs = [i for i in range(1, p.total_residue() + 1) if p.residue(i).is_ligand()]
    if not ligs:
        raise SystemExit(f"No ligand residue found in {pdb}")
    lig = ligs[-1]                     # works whatever chain/position the ligand has

    # binding-site residues: within 10 A of the ligand, SASA < 20 A^2 (single pass)
    sel_lig = rosetta.core.select.residue_selector.ResidueIndexSelector(str(lig))
    nb = rosetta.core.select.residue_selector.NeighborhoodResidueSelector(sel_lig, SHELL_CUTOFF, False)
    shell = []
    for r in rosetta.core.select.get_residues_from_subset(nb.apply(p)):
        s = rosetta.core.select.residue_selector.ResidueIndexSelector(str(r))
        m = rosetta.core.simple_metrics.metrics.SasaMetric()
        m.set_residue_selector(s)
        if m.calculate(p) < SHELL_SASA:
            shell.append(int(r))

    import json
    summary = {"design": stem, "ligand_residue": lig, "binding_site_residues": len(shell), "mutations": []}

    def finish(outcome, pose=None, energy_start=None):
        summary["outcome"] = outcome
        if pose is not None:
            summary["energy_change"] = round(sf(pose) - energy_start, 3) if energy_start is not None else 0.0
        with open(os.path.join(outdir, f"{stem}_hbrefine.json"), "w") as fh:
            json.dump(summary, fh, indent=1)
        return summary

    info = analyse(p, lig, shell)
    summary["unsat_before"] = [p.residue(lig).atom_name(a).strip() for a in info["unsat"]]
    e_start = sf(p)
    if not info["unsat"]:
        shutil.copy(pdb, os.path.join(outdir, f"{stem}_unrefined.pdb"))
        print("NO PROBLEMS: all buried ligand polar atoms are hydrogen bonded")
        summary["unsat_after"] = []
        return finish("already_satisfied", p, e_start)

    def repack(pose, keep):
        tp = rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
        tp.restrict_to_repacking()
        tp.temporarily_fix_everything()
        for r in shell:
            if r not in keep:
                tp.temporarily_set_pack_residue(r, True)
        rosetta.protocols.minimization_packing.PackRotamersMover(sf, tp).apply(pose)

    from pyrosetta.rosetta.core.scoring import ScoreType

    def minimize_site(pose):
        """Side chains (and backbone) of binding-site residues; ligand fixed."""
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        mm.set_jump(False)
        for r in shell:
            mm.set_chi(r, True)
            if min_backbone:
                mm.set_bb(r, True)
        mn = rosetta.protocols.minimization_packing.MinMover()
        mn.movemap(mm)
        mn.score_function(sf)
        mn.min_type("lbfgs_armijo_nonmonotone")
        mn.tolerance(0.001)
        mn.apply(pose)
        return sf(pose)

    def rep(pose):
        sf(pose)
        return pose.energies().total_energies()[ScoreType.fa_rep]

    nmut = 0
    for _ in range(len(info["unsat"]) + 1):
        info = analyse(p, lig, shell)
        if not info["unsat"]:
            break
        need = set()
        for role in info["roles"].values():
            if "acc" in role:
                need.update(DONOR_RES)
            if "don" in role:
                need.update(ACCEPTOR_RES)
        e0 = sf(p)
        rep0 = rep(p)
        base = None
        if minimize:                                  # baseline minimized the same way, so dE is not biased
            base = p.clone()
            e0_min = minimize_site(base)
        best = None
        diag = {"candidates": 0, "reducing": 0, "reducing_within_cap": 0, "best_reducing_dE": None,
                "best_reducing_mutation": None, "best_reducing_dE_fa_rep": None,
                "minimized": minimize, "best_reducing_dE_after_min": None}
        n_unsat = len(info["unsat"])
        best_min_pos = None
        for r in shell:
            if r in info["res_hb_lig"]:
                continue
            for aa in sorted(need):
                t = p.clone()
                mutate_residue(t, r, aa)
                repack(t, info["res_hb_lig"])
                e = sf(t)
                diag["candidates"] += 1
                if not (diagnose or minimize or e <= e0 + cap):
                    continue
                if len(analyse(t, lig, shell)["unsat"]) >= n_unsat:
                    continue                          # does not satisfy a buried ligand atom
                diag["reducing"] += 1
                mut = f"{p.residue(r).name1()}{r}{aa}"
                if diag["best_reducing_dE"] is None or e - e0 < diag["best_reducing_dE"]:
                    diag["best_reducing_dE"] = round(e - e0, 3)
                    diag["best_reducing_mutation"] = mut
                    diag["best_reducing_dE_fa_rep"] = round(rep(t) - rep0, 3)
                if minimize:
                    e_m = minimize_site(t)
                    if len(analyse(t, lig, shell)["unsat"]) >= n_unsat:
                        continue                      # H-bond lost during minimization
                    d_m = e_m - e0_min
                    if diag["best_reducing_dE_after_min"] is None or d_m < diag["best_reducing_dE_after_min"]:
                        diag["best_reducing_dE_after_min"] = round(d_m, 3)
                        diag["best_min_mutation"] = mut
                        best_min_pos = (r, aa)
                else:
                    e_m, d_m = e, e - e0
                if d_m <= cap:
                    diag["reducing_within_cap"] += 1
                    if best is None or e_m < best[0]:
                        best = (e_m, r, aa, t)
        if diagnose and minimize and best_min_pos is not None:
            # controls at the same position, with identical repacking and minimization:
            #   identity  - mutate the residue to itself (cost of the protocol alone)
            #   nonpolar  - size-matched nonpolar residue (same bulk, no hydrogen bond)
            r_c, aa_c = best_min_pos
            own = p.residue(r_c).name1()
            for tag, aa1 in (("self", own), ("np", NP_ANALOG.get(aa_c, "A"))):
                t = p.clone()
                mutate_residue(t, r_c, aa1)
                repack(t, info["res_hb_lig"])
                e_c = sf(t)
                e_cm = minimize_site(t)
                diag[f"control_{tag}_mutation"] = f"{own}{r_c}{aa1}"
                diag[f"control_{tag}_dE"] = round(e_c - e0, 3)
                diag[f"control_{tag}_dE_after_min"] = round(e_cm - e0_min, 3)
        if best is not None and minimize:
            best = (best[0] - e0_min + e0, best[1], best[2], best[3])   # report dE relative to the minimized baseline
        summary.setdefault("hbond_search", []).append(diag)
        print(f"H-bond step: {diag['candidates']} candidates, {diag['reducing']} reduce buried unsatisfied ligand atoms"
              + (f" (best dE {diag['best_reducing_dE']:+.2f} REU of which fa_rep {diag['best_reducing_dE_fa_rep']:+.2f}, "
                 f"{diag['best_reducing_mutation']})" if diag["reducing"] else "")
              + (f"; after minimization best dE {diag['best_reducing_dE_after_min']:+.2f} REU"
                 if diag["best_reducing_dE_after_min"] is not None else "")
              + f"; {diag['reducing_within_cap']} within the {cap} REU cap")
        if best is None:
            break
        print(f"H-bond mutation {p.residue(best[1]).name3()}{best[1]}{best[2]}: dE = {best[0] - e0:+.2f} REU")
        summary["mutations"].append(f"{p.residue(best[1]).name1()}{best[1]}{best[2]}")
        p = best[3]
        nmut += 1

    info = analyse(p, lig, shell)
    for r in useless_polar(p, lig, shell, info):
        e0 = sf(p)
        best = None
        for aa in NONPOLAR:
            t = p.clone()
            mutate_residue(t, r, aa)
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            mm.set_bb(r, True)
            mm.set_chi(r, True)
            mn = rosetta.protocols.minimization_packing.MinMover()
            mn.movemap(mm)
            mn.score_function(sf)
            mn.apply(t)
            e = sf(t)
            if e <= e0 + cap and (best is None or e < best[0]):
                best = (e, aa, t)
        if best is not None:
            print(f"buried polar mutation {p.residue(r).name3()}{r}{best[1]}: dE = {best[0] - e0:+.2f} REU")
            summary["mutations"].append(f"{p.residue(r).name1()}{r}{best[1]}")
            p = best[2]
            nmut += 1

    info = analyse(p, lig, shell)
    ubp = useless_polar(p, lig, shell, info)
    base = os.path.join(outdir, f"{stem}_{nmut}mutations")
    summary["unsat_after"] = [p.residue(lig).atom_name(a).strip() for a in info["unsat"]]
    summary["useless_polar_after"] = len(ubp)
    if not info["unsat"]:
        p.dump_pdb(base + "_fixed.pdb")
        return finish("fixed", p, e_start)
    if nmut > 0:
        p.dump_pdb(base + ("_bothproblems.pdb" if ubp else "_unsats.pdb"))
        return finish("partially_fixed", p, e_start)
    print("no solutions")
    return finish("no_solution", p, e_start)


if __name__ == "__main__":
    main()
