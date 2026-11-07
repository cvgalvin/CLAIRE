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
                tmp = dst.with_suffix(".tmp")
                tmp.write_bytes(data)
                tmp.rename(dst)
                return pid, src
            except Exception:
                continue
        return pid, "missing"

    workers = int(src_cfg.get("download_workers", 8))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        results = dict(ex.map(one, sorted(set(p.upper() for p in pdb_ids))))
    return root, results
