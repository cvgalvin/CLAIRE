"""Where structures and ligand information come from.

Replaces the requirement for a full local PDB-REDO mirror:

* Structures are fetched on demand, per PDB entry, into a shared cache (default ``~/.cache/claire``),
  from PDB-REDO, with optional fallback to the RCSB deposited model. A local PDB-REDO mirror is used
  directly if configured. Only the entries containing your fragments are downloaded, and the cache is
  reused across projects.
* The fragment -> PDB-ligand -> PDB-entry lookup can run against the RCSB search API (current), or
  against the original bundled 2022 database for exact reproduction of published runs.
"""
from __future__ import annotations

import json
import os
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .project import ClaireError, Project

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
PDBREDO_URL = "https://pdb-redo.eu/db/{id}/{id}_final.pdb"
RCSB_URL = "https://files.rcsb.org/download/{ID}.pdb"
UA = {"User-Agent": "claire-design (https://github.com/cvgalvin/CLAIRE)"}


def cache_dir() -> Path:
    d = Path(os.environ.get("CLAIRE_CACHE", Path.home() / ".cache" / "claire")) / "structures"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _get(url: str, data: bytes | None = None, timeout: int = 60, tries: int = 6) -> bytes:
    headers = dict(UA)
    if data is not None:
        headers["Content-Type"] = "application/json"
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=headers), timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404 or (400 <= e.code < 500 and e.code != 429):
                raise                      # client errors will not succeed on retry
            last = e
        except Exception as e:  # network hiccup
            last = e
        time.sleep(min(2 ** i, 30))       # 1, 2, 4, 8, 16, 30 s
    raise last  # type: ignore[misc]


# ---------------------------------------------------------------- ligand lookup
def _search(body: dict) -> list[str]:
    raw = _get(RCSB_SEARCH, json.dumps(body).encode())
    return [r["identifier"] for r in json.loads(raw).get("result_set", [])] if raw else []


def rcsb_ligands_with_substructure(smiles: str) -> list[str]:
    """Chemical component IDs (PDB ligand codes) that contain ``smiles`` as a substructure."""
    return _search({"query": {"type": "terminal", "service": "chemical",
                              "parameters": {"value": smiles, "type": "descriptor", "descriptor_type": "SMILES",
                                             "match_type": "sub-struct-graph-relaxed"}},
                    "return_type": "mol_definition", "request_options": {"return_all_hits": True}})


def rcsb_entries_with_ligand(comp_id: str) -> list[str]:
    return _search({"query": {"type": "terminal", "service": "text",
                              "parameters": {"attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                                             "operator": "exact_match", "value": comp_id}},
                    "return_type": "entry", "request_options": {"return_all_hits": True}})


def bundled_db(proj: Project) -> tuple[dict, dict]:
    d = proj.config["sources"].get("bundled_db")
    if not d:
        raise ClaireError("sources.ligand_search is 'bundled' but sources.bundled_db is not set. Unzip db.zip from the "
                          "original repository and point sources.bundled_db at the db/ folder.")
    d = proj.path(d)
    smiles = json.loads((d / "PDB_small_molecule_SMILES_Inchi.json").read_text())
    lig2pdb = {}
    for line in (d / "pdb_small_molecules_ligandexpo.txt").read_text().splitlines():
        parts = [p for p in line.replace("\t", " ").split(" ") if p]
        if parts:
            lig2pdb[parts[0]] = parts[1:]
    return smiles, lig2pdb


# ---------------------------------------------------------------- local Chemical Component Dictionary
CCD_URL = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
_PREF = [("SMILES_CANONICAL", "OpenEye"), ("SMILES_CANONICAL", "CACTVS"), ("SMILES", "OpenEye"), ("SMILES", "CACTVS"),
         ("SMILES", "ACDLabs")]


def _parse_ccd(gz_path: Path) -> dict:
    """{comp_id: [(type, program, smiles), ...]} from components.cif.gz (descriptor loop only)."""
    import gzip
    import shlex
    out, comp, in_loop, pending = {}, None, False, None
    with gzip.open(gz_path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("data_"):
                comp, in_loop, pending = line[5:].strip(), False, None
                continue
            if line.startswith("_pdbx_chem_comp_descriptor."):
                in_loop = True
                continue
            if not in_loop:
                continue
            if line.startswith(("loop_", "#", "_")) and pending is None:
                in_loop = line.startswith("_pdbx_chem_comp_descriptor.")
                continue
            if pending is not None:                       # multi-line ';' value
                if line.startswith(";"):
                    if pending[3]:
                        out.setdefault(comp, []).append((pending[0], pending[1], pending[3].strip()))
                        pending = None
                    else:
                        pending[3] = line[1:].strip()
                else:
                    pending[3] += line.strip()
                continue
            try:
                parts = shlex.split(line)
            except ValueError:
                continue
            if len(parts) >= 5 and parts[0] == comp:
                out.setdefault(comp, []).append((parts[1], parts[2], parts[4]))
            elif len(parts) == 4 and parts[0] == comp:
                pending = [parts[1], parts[2], parts[3], ""]
    return out


def ccd_smiles(max_age_days: float = 30.0) -> dict:
    """{comp_id: SMILES} for every chemical component in the current CCD (downloaded and cached)."""
    root = Path(os.environ.get("CLAIRE_CACHE", Path.home() / ".cache" / "claire"))
    root.mkdir(parents=True, exist_ok=True)
    js = root / "ccd_smiles.json"
    if js.exists() and (time.time() - js.stat().st_mtime) < max_age_days * 86400:
        return json.loads(js.read_text())
    gz = root / "components.cif.gz"
    print("[fragments] downloading the Chemical Component Dictionary (~120 MB, cached for later runs)")
    tmp = gz.with_suffix(f".{os.getpid()}.tmp")
    tmp.write_bytes(_get(CCD_URL, timeout=600))
    os.replace(tmp, gz)
    desc = _parse_ccd(gz)
    smiles = {}
    for comp, rows in desc.items():
        for t, prog in _PREF:
            hit = next((d for tt, pp, d in rows if tt == t and pp.startswith(prog)), None)
            if hit:
                smiles[comp] = hit
                break
    tmpj = js.with_suffix(f".{os.getpid()}.tmp")
    tmpj.write_text(json.dumps(smiles))
    os.replace(tmpj, js)
    return smiles


# ---------------------------------------------------------------- structures
def _layout(root: Path, pdbid: str) -> Path:
    i = pdbid.lower()
    return root / i[1:3] / i / f"{i}_final.pdb"


def fetch_structures(proj: Project, pdb_ids: list[str]) -> tuple[Path, dict]:
    """Ensure every entry is available under the cache, in the PDB-REDO directory layout.

    Returns (cache root, {pdb_id: source}) where source is 'local', 'cache', 'pdb-redo', 'rcsb' or 'missing'.
    """
    src_cfg = proj.config["sources"]
    root = cache_dir()
    mirror = proj.path(src_cfg.get("local_pdb_redo")) if src_cfg.get("local_pdb_redo") else None
    prefer = src_cfg.get("structures", "pdb-redo")
    fallback = bool(src_cfg.get("fallback_to_rcsb", True))

    def one(pid: str) -> tuple[str, str]:
        dst = _layout(root, pid)
        if dst.exists() and dst.stat().st_size > 0:
            return pid, "cache"
        dst.parent.mkdir(parents=True, exist_ok=True)
        if mirror is not None:
            m = _layout(mirror, pid)
            if m.exists():
                try:
                    os.symlink(m, dst)
                except FileExistsError:
                    pass
                return pid, "local"
        order = ["pdb-redo", "rcsb"] if prefer == "pdb-redo" else ["rcsb"]
        if prefer == "pdb-redo" and not fallback:
            order = ["pdb-redo"]
        for src in order:
            url = PDBREDO_URL.format(id=pid.lower()) if src == "pdb-redo" else RCSB_URL.format(ID=pid.upper())
            try:
                data = _get(url, timeout=90)
                if data.lstrip()[:6] in (b"<!DOCT", b"<html>") or len(data) < 200:
                    continue
                tmp = dst.with_suffix(f".{os.getpid()}.{threading.get_ident()}.tmp")   # unique: safe with shared caches
                tmp.write_bytes(data)
                os.replace(tmp, dst)
                return pid, src
            except Exception:
                continue
        return pid, "missing"

    workers = int(src_cfg.get("download_workers", 8))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        results = dict(ex.map(one, sorted(set(p.upper() for p in pdb_ids))))
    return root, results
