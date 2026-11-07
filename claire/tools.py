"""Locate external programs (Rosetta, ProteinMPNN, ColabFold, Open Babel) and packaged resources."""
from __future__ import annotations

import glob
import os
import shutil
import sys
from importlib import resources
from pathlib import Path

from .project import ClaireError, Project


def xml(name: str) -> str:
    """Absolute path of a packaged RosettaScripts XML (e.g. 'fd_mpnn.xml' or 'extra/scaffold_metrics.xml')."""
    return str(resources.files("claire.data").joinpath("xml", name))


def legacy(name: str) -> str:
    """Absolute path of a vendored original CLAIRE script."""
    return str(resources.files("claire.legacy").joinpath(name))


def python(proj: Project) -> str:
    return proj.config["tools"].get("python") or sys.executable


def rosetta_bin_dir(proj: Project) -> Path:
    b = proj.config["tools"].get("rosetta_bin") or os.environ.get("ROSETTA_BIN")
    if not b:
        raise ClaireError("Rosetta not configured: set tools.rosetta_bin in claire.yaml "
                          "(…/rosetta/main/source/bin) or the ROSETTA_BIN environment variable.")
    p = proj.path(b)
    if not p.is_dir():
        raise ClaireError(f"tools.rosetta_bin does not exist: {p}")
    return p


def rosetta(proj: Project, app: str) -> str:
    """Full path to a Rosetta application, e.g. rosetta(proj, 'match')."""
    bindir = rosetta_bin_dir(proj)
    suffix = proj.config["tools"].get("rosetta_suffix", "auto")
    if suffix and suffix != "auto":
        exe = bindir / f"{app}{suffix}"
        if exe.exists():
            return str(exe)
        raise ClaireError(f"Rosetta application not found: {exe}")
    hits = sorted(glob.glob(str(bindir / f"{app}.*release"))) or sorted(glob.glob(str(bindir / f"{app}.*")))
    if not hits:
        raise ClaireError(f"No Rosetta '{app}' binary in {bindir}. Set tools.rosetta_suffix explicitly.")
    # prefer default/static builds over mpi
    hits.sort(key=lambda h: ("mpi" in h, len(h)))
    return hits[0]


def molfile_to_params(proj: Project) -> str:
    p = proj.config["tools"].get("molfile_to_params")
    if p:
        return str(proj.path(p))
    guess = rosetta_bin_dir(proj).parent / "scripts" / "python" / "public" / "molfile_to_params.py"
    if guess.exists():
        return str(guess)
    raise ClaireError("molfile_to_params.py not found; set tools.molfile_to_params in claire.yaml.")


def mpnn(proj: Project) -> tuple[str, Path]:
    """Locate the MPNN checkout. Returns ('ligandmpnn', dir) for a LigandMPNN repository (run.py) or
    ('proteinmpnn', dir) for the original ProteinMPNN repository (protein_mpnn_run.py)."""
    d = proj.config["tools"].get("proteinmpnn_dir")
    if not d:
        raise ClaireError("Set tools.proteinmpnn_dir in claire.yaml to a ProteinMPNN or LigandMPNN checkout.")
    d = proj.path(d)
    if (d / "run.py").exists() and (d / "model_params").exists() or (d / "run.py").exists() and (d / "get_model_params.sh").exists():
        return "ligandmpnn", d
    for cand in (d, d / "vanilla_proteinmpnn"):
        if (cand / "protein_mpnn_run.py").exists():
            return "proteinmpnn", cand
    raise ClaireError(f"{d} is neither a LigandMPNN checkout (run.py) nor a ProteinMPNN checkout (protein_mpnn_run.py).")


LIGANDMPNN_WEIGHTS = {"protein_mpnn": ("checkpoint_protein_mpnn", "proteinmpnn_v_48_020.pt"),
                      "ligand_mpnn": ("checkpoint_ligand_mpnn", "ligandmpnn_v_32_010_25.pt")}


def ligandmpnn_checkpoint(proj: Project, d: Path, model_type: str) -> tuple[str, Path]:
    if model_type not in LIGANDMPNN_WEIGHTS:
        raise ClaireError(f"mpnn.model_type must be one of {sorted(LIGANDMPNN_WEIGHTS)}")
    flag, fname = LIGANDMPNN_WEIGHTS[model_type]
    ck = d / "model_params" / fname
    if not ck.exists():
        raise ClaireError(f"LigandMPNN weights not found ({ck}). Run `bash get_model_params.sh \"./model_params\"` in {d}.")
    return flag, ck


def executable(proj: Project, key: str) -> str:
    exe = proj.config["tools"].get(key) or key
    found = shutil.which(exe) or (exe if os.path.exists(exe) else None)
    if not found:
        raise ClaireError(f"'{exe}' not found on PATH; install it or set tools.{key} in claire.yaml.")
    return found
