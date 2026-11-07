"""`claire check`: report which dependencies are available, and which steps they are needed for."""
from __future__ import annotations

import importlib
import shutil

from . import tools
from .project import ClaireError, Project

PY = [("numpy", "all"), ("scipy", "contacts, cluster"), ("pandas", "contacts, cluster"), ("matplotlib", "cluster"),
      ("yaml", "all"), ("rdkit", "target, fragments, contacts"), ("prody", "contacts"), ("Bio", "contacts"),
      ("pyrosetta", "cluster (H-bond screen), motifs, enzdes, hbrefine")]


def _ok(flag: bool) -> str:
    return "ok     " if flag else "MISSING"


def _mpnn_desc(proj: Project) -> str:
    kind, d = tools.mpnn(proj)
    model = proj.config["mpnn"].get("model_type", "protein_mpnn")
    if kind == "ligandmpnn":
        _, ck = tools.ligandmpnn_checkpoint(proj, d, model)
        return f"{d} (LigandMPNN, model {model}, weights {ck.name})"
    return f"{d} (ProteinMPNN)"


def run_check(proj: Project) -> int:
    print("Python packages")
    missing = 0
    for mod, used in PY:
        try:
            importlib.import_module(mod)
            print(f"  {_ok(True)} {mod:<11} ({used})")
        except ImportError:
            missing += 1
            print(f"  {_ok(False)} {mod:<11} ({used})")
    print("External tools")
    checks = [
        ("Rosetta match", lambda: tools.rosetta(proj, "match"), "match"),
        ("Rosetta rosetta_scripts", lambda: tools.rosetta(proj, "rosetta_scripts"), "burial, score, fastdesign"),
        ("Rosetta enzyme_design", lambda: tools.rosetta(proj, "enzyme_design"), "enzdes"),
        ("molfile_to_params.py", lambda: tools.molfile_to_params(proj), "target"),
        ("ProteinMPNN/LigandMPNN", lambda: _mpnn_desc(proj), "mpnn"),
        ("colabfold_batch", lambda: tools.executable(proj, "colabfold_batch"), "predict"),
        ("obabel (optional)", lambda: tools.executable(proj, "obabel"), "fragments"),
    ]
    import subprocess
    for name, fn, used in checks:
        try:
            where = fn()
            if name.startswith("Rosetta"):            # make sure the binary actually launches
                r = subprocess.run([where, "-help"], capture_output=True, text=True, timeout=60)
                out = (r.stdout + r.stderr)
                if "Library not loaded" in out or "cannot open shared object" in out:
                    missing += 1
                    lib = next((l.strip() for l in out.splitlines() if "Library not loaded" in l or "shared object" in l), "")
                    print(f"  {_ok(False)} {name:<24} {where}\n             cannot start: {lib}\n"
                          "             The Rosetta folder was probably moved after compiling; symlink the original path "
                          "(shown above) to the current Rosetta folder, or recompile.")
                    continue
            print(f"  {_ok(True)} {name:<24} {where}")
        except subprocess.TimeoutExpired:
            print(f"  {_ok(True)} {name:<24} {where} (started; -help timed out)")
        except ClaireError as e:
            missing += "optional" not in name
            print(f"  {_ok(False)} {name:<24} ({used}) {e}")
    print("\nAll dependencies found." if not missing else
          f"\n{missing} missing; steps that need them will stop with an error. See README 'Installation'.")
    return 0
