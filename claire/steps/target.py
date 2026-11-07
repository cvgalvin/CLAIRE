"""Step 1: Rosetta params file and model for the target ligand (optionally with RDKit conformers)."""
from __future__ import annotations

import shutil
import subprocess

from .. import tools
from ..project import ClaireError, Project


def run(proj: Project, args) -> None:
    cfg = proj.config["ligand"]
    sdf = proj.path(cfg["sdf"])
    if not sdf.exists():
        raise ClaireError(f"Ligand model not found: {sdf} (set ligand.sdf in claire.yaml).")
    name = proj.lig
    if len(name) > 3:
        raise ClaireError("ligand.name must be at most 3 characters (Rosetta residue name).")
    out = proj.stage("target")
    local = out / f"{name}.sdf"
    shutil.copy(sdf, local)
    extra = []
    n = int(cfg.get("conformers", 1))
    if n > 1:
        print(f"[target] generating up to {n} conformers with RDKit")
        subprocess.run([tools.python(proj), tools.legacy("gencon.py"), local.name, str(n), "1000",
                        str(cfg.get("conformer_rms_prune", 0.1))], cwd=out, check=True)
        conf = out / "1.sdf"
        if not conf.exists():
            raise ClaireError("Conformer generation produced no output (see messages above).")
        conf.rename(out / f"{name}_conformers.sdf")
        local = out / f"{name}_conformers.sdf"
        extra = ["--conformers-in-one-file"]
    cmd = [tools.python(proj), tools.molfile_to_params(proj), "-n", name, "-p", name, *extra, local.name]
    print("[target]", " ".join(cmd))
    r = subprocess.run(cmd, cwd=out, capture_output=True, text=True)
    (out / "molfile_to_params.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0 or not (out / f"{name}.params").exists():
        raise ClaireError(f"molfile_to_params failed; see {out / 'molfile_to_params.log'}")
    proj.write_manifest("target", {"params": f"{name}.params", "pdb": f"{name}_0001.pdb",
                                   "conformers": n, "n_outputs": 1})
    print(f"[target] wrote {out / (name + '.params')} and {out / (name + '_0001.pdb')}")
    print("Next: define fragments in claire.yaml (see `claire fragments --list-atoms`), then run `claire fragments`.")
