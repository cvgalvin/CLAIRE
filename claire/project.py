"""Project layout, configuration and per-stage manifests.

A CLAIRE project is a directory containing ``claire.yaml`` plus one numbered folder per workflow
stage. Each stage writes a ``manifest.json`` describing its outputs, which the next stage reads,
so no file paths have to be passed between steps by hand.
"""
from __future__ import annotations

import copy
import json
import os
import shutil
from importlib import resources
from pathlib import Path

import yaml

CONFIG_NAME = "claire.yaml"

# stage key -> folder name (order = workflow order)
STAGES = {
    "target": "01_target",
    "fragments": "02_fragments",
    "contacts": "03_contacts",
    "cluster": "04_clusters",
    "motifs": "05_motifs",
    "match": "06_matches",
    "burial": "07_buried",
    "enzdes": "08_enzdes",
    "hbrefine": "09_hbrefine",
    "score": "10_scored",
    "mpnn": "11_mpnn",
    "fastdesign": "12_fastdesign",
    "predict": "13_colabfold",
    "final": "14_final",
}


class ClaireError(RuntimeError):
    """User-facing error with an actionable message."""


def template_config() -> dict:
    text = resources.files("claire.data").joinpath("config_template.yaml").read_text()
    return yaml.safe_load(text)


def _merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


class Project:
    def __init__(self, root: str | Path):
        self.root = Path(root).resolve()
        cfg_path = self.root / CONFIG_NAME
        if not cfg_path.exists():
            raise ClaireError(f"No {CONFIG_NAME} in {self.root}. Run `claire init {self.root}` first, "
                              "or run commands from inside a project directory (or pass --project).")
        try:
            user = yaml.safe_load(cfg_path.read_text()) or {}
        except yaml.YAMLError as e:
            mark = getattr(e, "problem_mark", None)
            where = f" at line {mark.line + 1}, column {mark.column + 1}" if mark else ""
            line = ""
            if mark:
                lines = cfg_path.read_text().splitlines()
                prev = [i for i in range(mark.line - 1, -1, -1)
                        if lines[i].strip() and not lines[i].lstrip().startswith("#")][:3]
                ctx = sorted(prev) + ([mark.line] if mark.line < len(lines) else [])
                line = "".join(f"\n    {i + 1:>4} | {lines[i]}" for i in ctx) + "\n  (the mistake is usually on or just above the last line shown)"
            raise ClaireError(f"{cfg_path.name} is not valid YAML{where}:{line}\n"
                              f"  ({getattr(e, 'problem', e)}). Each setting is `key: value` on its own line, indented "
                              "under its section; {key: value} braces only work as the whole value of a section, "
                              "e.g. `enzdes: {filter: false}`.")
        self.config = _merge(template_config(), user)

    # ---------- discovery ----------
    @classmethod
    def find(cls, start: str | Path | None = None) -> "Project":
        p = Path(start or os.getcwd()).resolve()
        for d in [p, *p.parents]:
            if (d / CONFIG_NAME).exists():
                return cls(d)
        raise ClaireError(f"Could not find {CONFIG_NAME} in {p} or any parent directory. "
                          "Run `claire init <dir>` to create a project.")

    # ---------- paths ----------
    def path(self, rel: str | None) -> Path | None:
        if rel is None:
            return None
        rel = os.path.expandvars(os.path.expanduser(str(rel)))
        p = Path(rel)
        return p if p.is_absolute() else (self.root / p)

    def stage(self, key: str, create: bool = True) -> Path:
        d = self.root / STAGES[key]
        if create:
            d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def lig(self) -> str:
        return str(self.config["ligand"]["name"])

    def params(self) -> Path:
        p = self.stage("target", create=False) / f"{self.lig}.params"
        if not p.exists():
            raise ClaireError("Ligand params not found. Run `claire target` first.")
        return p

    def ligand_pdb(self) -> Path:
        p = self.stage("target", create=False) / f"{self.lig}_0001.pdb"
        if not p.exists():
            raise ClaireError("Ligand model not found. Run `claire target` first.")
        return p

    # ---------- fragments ----------
    def fragments(self) -> list[dict]:
        """Normalized fragment definitions: [{'name': str, 'atoms': [str, ...]}, ...].

        Accepted forms in claire.yaml:
          fragments: [C1, C2, C3]                     -> one fragment
          fragments: [[C1, C2, C3], [C17, C20, O20]]  -> two fragments
          fragments: [{name: ring, atoms: [C1, C2, C3]}, ...]
        """
        raw = self.config.get("fragments") or []
        if raw == "auto":
            raw = [{"name": f["name"], "atoms": f["atoms"]} for f in self.auto_fragments()["fragments"]]
        if raw and all(isinstance(x, str) for x in raw):
            raw = [raw]
        out = []
        for i, f in enumerate(raw, 1):
            if isinstance(f, dict):
                atoms, name = f.get("atoms"), f.get("name", f"Fragment_{i}")
            elif isinstance(f, (list, tuple)):
                atoms, name = list(f), f"Fragment_{i}"
            else:
                atoms = None
            if not atoms or not all(isinstance(a, str) for a in atoms):
                raise ClaireError(
                    f"fragments entry {i} is not understood: {f!r}. Use, for example:\n"
                    "fragments:\n  - name: ring\n    atoms: [C1, C2, C3, C4, C5, C10]\n"
                    "or the shorthand  fragments: [[C1, C2, C3, C4, C5, C10], [C17, C20, O20, C21]]")
            out.append({"name": str(name), "atoms": [a.strip() for a in atoms]})
        return out

    def auto_fragments(self) -> dict:
        """Automatic fragment definitions (cached in 02_fragments/auto_fragments.yaml)."""
        f = self.root / STAGES["fragments"] / "auto_fragments.yaml"
        if f.exists():
            return yaml.safe_load(f.read_text())
        from .steps.autofrag import define
        d = define(self)
        d.pop("roles", None)
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(yaml.safe_dump(d, sort_keys=False))
        return d

    def hbond_fragments(self) -> list[str]:
        hb = self.config.get("hbond_fragments") or []
        if hb == "auto":
            return list(self.auto_fragments()["hbond_fragments"])
        return list(hb)

    # ---------- manifests ----------
    def write_manifest(self, key: str, data: dict) -> None:
        data = dict(data)
        data["stage"] = key
        (self.stage(key) / "manifest.json").write_text(json.dumps(data, indent=1, default=str))

    def manifest(self, key: str, required: bool = True) -> dict:
        p = self.root / STAGES[key] / "manifest.json"
        if not p.exists():
            if required:
                raise ClaireError(f"Stage '{key}' has not been completed ({p} missing). Run `claire {key}` first.")
            return {}
        return json.loads(p.read_text())

    def status(self) -> list[tuple[str, str, str]]:
        rows = []
        for key, folder in STAGES.items():
            m = self.manifest(key, required=False)
            if m:
                n = m.get("n_outputs", len(m.get("outputs", [])))
                rows.append((key, folder, f"done ({n} outputs)"))
            elif (self.root / folder / "jobs").exists():
                rows.append((key, folder, "submitted; run with --collect when jobs finish"))
            else:
                rows.append((key, folder, "-"))
        return rows


def init_project(root: str | Path, ligand_sdf: str | None, name: str | None, scaffolds: str | None) -> Path:
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    (root / "inputs" / "scaffolds").mkdir(parents=True, exist_ok=True)
    cfg_path = root / CONFIG_NAME
    if cfg_path.exists():
        raise ClaireError(f"{cfg_path} already exists.")
    text = resources.files("claire.data").joinpath("config_template.yaml").read_text()
    if name:
        text = text.replace("name: LIG ", f"name: {name.upper()} ", 1)
    if ligand_sdf:
        dst = root / "inputs" / "ligand.sdf"
        shutil.copy(ligand_sdf, dst)
    if scaffolds:
        src = Path(scaffolds)
        n = 0
        for f in src.iterdir():
            if f.suffix.lower() in (".pdb", ".pos"):
                shutil.copy(f, root / "inputs" / "scaffolds" / f.name)
                n += f.suffix.lower() == ".pdb"
        if n == 0:
            raise ClaireError(f"No .pdb files found in {src}")
    cfg_path.write_text(text)
    return cfg_path
