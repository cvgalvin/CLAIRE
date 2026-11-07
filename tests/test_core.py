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
    b = motifs._hbond_block("SER", "STR", ["O3", "C3", "C2"], False, "y")
    assert "distanceAB:    1.90   0.50" in b[4] and "angle_A:  120.00  10.00" in b[5]


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
