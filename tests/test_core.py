"""Tests that run without Rosetta, PyRosetta, ProteinMPNN or ColabFold."""
import argparse
import numpy as np
import pytest

from claire.project import Project, init_project, ClaireError
from claire.steps import motifs, design, predict
from claire.executor import Task, run_tasks


@pytest.fixture
def proj(tmp_path):
    init_project(tmp_path / "p", None, "STR", None)
    cfg = (tmp_path / "p" / "claire.yaml")
    cfg.write_text(cfg.read_text().replace("polar_atoms: []\n",
        "polar_atoms:\n  - atoms: [O3, C3, C2]\n  - atoms: [O20, C20, C21]\n", 1))
    return Project(tmp_path / "p")


def test_config_defaults(proj):
    assert proj.lig == "STR"
    assert proj.config["burial"]["dsasa_min"] == 0.7
    assert proj.config["final"]["plddt_min"] == 88.0


def test_missing_stage_is_actionable(proj):
    with pytest.raises(ClaireError, match="claire target"):
        proj.params()


def test_polar_motifs_follow_original_rules(proj):
    import random
    random.seed(1)
    m = motifs.polar_motifs(proj)
    assert m and all(len(x) == 2 for x in m)
    assert len({str(x) for x in m}) == len(m)                       # unique
    for x in m:
        res = [[l for l in b if "ATOM_MAP: 2 residue3" in l][0].split()[-1] for b in x]
        assert res.count("GLY") <= 2
        assert sum(r in ("LYS", "ARG") for r in res) <= 1


def test_sp2_acceptor_block_values():
    motifs.set_overrides({}, {})
    motifs.GEOMETRY["set"] = "original"
    b = motifs._hbond_block("SER", "STR", ["O3", "C3", "C2"], False, "y")
    assert "distanceAB:    1.90   0.50" in b[4] and "angle_A:  120.00  10.00" in b[5]
    motifs.GEOMETRY["set"] = "revised"
    b = motifs._hbond_block("SER", "STR", ["O3", "C3", "C2"], False, "y")
    assert "distanceAB:    2.00   0.30" in b[4] and "angle_A:  135.00  20.00" in b[5]
    assert b[7].rstrip().endswith(" 2") and b[8].rstrip().endswith(" 9")   # torsion_A constrained, torsion_B free


def test_donor_requires_hydrogen_first(proj):
    proj.config["polar_atoms"] = [{"atoms": ["O1", "C3", "C2"], "donor": True}]
    with pytest.raises(ClaireError, match="polar hydrogen"):
        motifs.polar_motifs(proj)


def test_motif_residue_parsing():
    assert design.motif_residues("UM_12_Y14T54T98M112_x__1np__hybrid_3.pdb") == [14, 54, 98, 112]


def test_ca_rmsd_superposition_invariant():
    x = np.random.RandomState(0).rand(30, 3) * 10
    rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert predict.ca_rmsd(x, x @ rot.T + 3) < 1e-8


def test_local_executor(proj, tmp_path):
    tasks = [Task(f"echo {i} > out{i}", str(tmp_path), f"t{i}") for i in range(3)]
    assert run_tasks(proj, "match", tasks)
    assert (tmp_path / "out2").read_text().strip() == "2"


def test_slurm_script_written(proj, tmp_path):
    proj.config["compute"]["backend"] = "slurm"
    assert not run_tasks(proj, "match", [Task("true", str(tmp_path), "t")], label="demo")
    assert "#SBATCH --array=1-1" in (proj.root / "06_matches" / "jobs" / "demo.sh").read_text()


def test_constraint_overrides_and_fast_preset(proj):
    motifs.set_overrides({}, {})
    default = motifs._hbond_block("SER", "STR", ["O3", "C3", "C2"], False, "y")
    fp = proj.config["motifs"]["fast_preset"]
    motifs.set_overrides(fp["hbond"], fp["nonpolar"])
    fast = motifs._hbond_block("SER", "STR", ["O3", "C3", "C2"], False, "y")
    motifs.set_overrides({}, {})
    assert default != fast
    assert all(l.rstrip().endswith(" 0") for l in fast if "CONSTRAINT" in l)          # ideal values only
    assert any("360.00 100.00" in l for l in fast if "torsion_B" in l)                 # free torsion stays free
    with pytest.raises(ClaireError):
        motifs.set_overrides({"not_a_key": 1}, {})


def _write_ligand(proj, smiles, name):
    """Write a params-style file and model with simple atom names for a small test ligand."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=1)
    names, cnt = [], {}
    for a in m.GetAtoms():
        e = a.GetSymbol()
        cnt[e] = cnt.get(e, 0) + 1
        names.append(f"{e}{cnt[e]}")
    d = proj.root / "01_target"
    d.mkdir(exist_ok=True)
    pos = m.GetConformer().GetPositions()
    (d / f"{name}_0001.pdb").write_text("".join(
        f"HETATM{i + 1:5d} {n:<4s} {name} X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a.GetSymbol():>2s}\n"
        for i, (n, a, (x, y, z)) in enumerate(zip(names, m.GetAtoms(), pos))))
    bonds = "".join(f"BOND_TYPE {names[b.GetBeginAtomIdx()]} {names[b.GetEndAtomIdx()]} 1\n" for b in m.GetBonds())
    (d / f"{name}.params").write_text(f"NAME {name}\n" + "".join(f"ATOM  {n} X X 0.0\n" for n in names) + bonds)


def test_auto_fragments_phenol_and_ketone(proj):
    from claire.steps import autofrag
    proj.config["ligand"]["name"] = "LIG"
    # phenol ring + ketone + a separate benzyl ring; params bond orders all single
    _write_ligand(proj, "Oc1ccc(CC(C)=O)cc1", "LIG")
    d = autofrag.define(proj)
    frags = {f["name"]: set(f["atoms"]) for f in d["fragments"]}
    phenol = [f for n, f in frags.items() if n.startswith("hb_") and "O1" in f][0]
    assert {"C1", "C2", "C3", "C4", "C8", "C9"} <= phenol             # whole ring with the OH (orders perceived)
    assert any(n.startswith("hb_") and "O2" in f and "C6" in f for n, f in frags.items())   # ketone found despite all-single params
    assert not any(n.startswith("aromatic_") for n in frags)          # ring already covered by the phenol fragment
    proj.config["ligand"]["name"] = "LIH"
    _write_ligand(proj, "Oc1ccc(CCCc2ccccc2)cc1", "LIH")             # second ring not covered by any polar group
    d2 = autofrag.define(proj)
    assert sum(f["name"].startswith("aromatic_") for f in d2["fragments"]) == 1
    assert any(s["donor"] and s["atoms"][0].startswith("H") for s in d["polar_atoms_suggested"])


def test_statistical_fit():
    from claire.steps import hbstats
    rng = np.random.RandomState(0)
    ideal, tol, n = hbstats._fit(rng.normal(150, 5, 200), "angle")
    assert abs(ideal - 150) < 2 and 10 <= tol <= 35
    assert hbstats._fit(rng.uniform(-180, 180, 500), "torsion")[1] == 360      # no preference -> free
    t = hbstats._fit(np.r_[rng.normal(178, 8, 100), rng.normal(-178, 8, 100)], "torsion")
    assert abs(abs(t[0]) - 180) < 5                                              # circular mean across +/-180


def test_fast_clustering_matches_original(tmp_path):
    import filecmp
    import os
    import subprocess
    import sys
    from claire import tools  # noqa: F401
    from importlib import resources
    rng = np.random.RandomState(3)
    names = ["N", "CA", "C", "O", "CB", "CG", "OG1", "CG2"]
    modes = [rng.normal(0, 3, (8, 3)) for _ in range(5)]
    lines = []
    for r in range(1, 301):
        xyz = modes[rng.randint(5)] + rng.normal(0, 0.25, (8, 3))
        keep = range(6) if r % 17 == 0 else range(8)                       # some truncated side chains
        for k in keep:
            l = f"ATOM  {1:5d} {names[k]:<4s} THR A{0:4d}    {xyz[k,0]:8.3f}{xyz[k,1]:8.3f}{xyz[k,2]:8.3f}  1.00  0.00           {names[k][0]}\n"
            lines.append(l[0:22] + (" " + str(r)).ljust(7) + l[29:])
    lig = tmp_path / "lig.pdb"
    lig.write_text("HETATM    1 C1   LIG X   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n")
    outs = {}
    for mode in ("orig", "fast"):
        d = tmp_path / mode
        d.mkdir()
        (d / "X_THR_contact_fuzzball.pdb").write_text("".join(lines))
        cmd = ([sys.executable, str(resources.files("claire.legacy").joinpath("spatial_cluster.py"))] if mode == "orig"
               else [sys.executable, "-m", "claire.cluster_fast"])
        subprocess.run(cmd + ["X_THR_contact_fuzzball.pdb", str(lig)], cwd=d, check=True, capture_output=True)
        outs[mode] = d / "X_THR_contact_fuzzball_clustered"
    fa, fb = sorted(os.listdir(outs["orig"])), sorted(os.listdir(outs["fast"]))
    assert fa == fb and len(fa) > 5
    assert all(filecmp.cmp(outs["orig"] / x, outs["fast"] / x, shallow=False) for x in fa)


def test_drop_impossible_contacts(tmp_path):
    from claire.steps.contacts import drop_impossible_contacts
    lig = ("HETATM    1  C1  LIG X   1       0.000   0.000   0.000  1.00  0.00           C\n"
           "HETATM    2  O1  LIG X   1       1.200   0.000   0.000  1.00  0.00           O\n")
    hbond = "ATOM      1  OG  SER A   1       3.800   0.000   0.000  1.00  0.00           O\n"   # O...O 2.6 A: keep
    clash = "ATOM      2  CB  VAL B   2       0.000   2.500   0.000  1.00  0.00           C\n"   # C...C 2.5 A: remove
    close = "ATOM      3  OG  SER C   3       3.500   0.000   0.000  1.00  0.00           O\n"   # O...O 2.3 A: remove
    pack = "ATOM      4  CB  LEU D   4       0.000   3.600   0.000  1.00  0.00           C\n"    # C...C 3.6 A: keep
    f = tmp_path / "x.pdb"
    f.write_text(hbond + clash + close + pack + lig)
    kept, removed = drop_impossible_contacts(f, {"C1", "O1"}, "LIG")
    assert (kept, removed) == (2, 2)
    text = f.read_text()
    assert "SER A" in text and "LEU D" in text and "VAL B" not in text and "SER C" not in text and "LIG X" in text
