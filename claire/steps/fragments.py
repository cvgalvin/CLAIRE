"""Step 2: build fragment models from ligand atom names and find PDB entries containing each fragment.

Output (02_fragments/): Fragment_<n>.pdb/.mol and PDB_search_results.json, in the same format as the
original Fragment_Search.py so the vendored extraction script can consume it unchanged.
"""
from __future__ import annotations

import json
import subprocess
import urllib.request
from pathlib import Path

from .. import sources, tools
from ..project import ClaireError, Project


COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39, "B": 0.84, "SI": 1.11}


def _ligand_atoms(proj: Project):
    """(name, element, xyz) for each atom of the ligand model."""
    out = []
    for l in proj.ligand_pdb().read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            name = l[12:16].strip()
            el = (l[76:78].strip() or "".join(c for c in name if c.isalpha())[:1]).upper()
            out.append((name, el, (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return out


def _bonds(atoms) -> dict:
    """Bonded neighbours from interatomic distances (covalent radii + 0.4 A)."""
    import math
    nb = {a[0]: [] for a in atoms}
    for i, (n1, e1, x1) in enumerate(atoms):
        for n2, e2, x2 in atoms[i + 1:]:
            if math.dist(x1, x2) <= COV.get(e1, 0.8) + COV.get(e2, 0.8) + 0.4:
                nb[n1].append(n2)
                nb[n2].append(n1)
    return nb


def list_atoms(proj: Project) -> None:
    atoms = _ligand_atoms(proj)
    nb = _bonds(atoms)
    print(f"Atoms in {proj.ligand_pdb().name} and their bonded neighbours (fragment atoms must form one connected piece):")
    for name, el, _ in atoms:
        heavy = [n for n in nb[name] if not n.upper().startswith("H")]
        hyd = [n for n in nb[name] if n.upper().startswith("H")]
        if el != "H":
            print(f"  {name:<5} {el:<2}  bonded to {', '.join(heavy) or '-'}" + (f"   (+H: {', '.join(hyd)})" if hyd else ""))


def _check_connected(proj: Project, fname: str, atoms: list[str]) -> None:
    nb = _bonds(_ligand_atoms(proj))
    sel = set(atoms)
    seen, parts = set(), []
    for a in atoms:
        if a in seen:
            continue
        comp, stack = [], [a]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.append(x)
            stack += [y for y in nb.get(x, []) if y in sel and y not in seen]
        parts.append(comp)
    if len(parts) > 1:
        desc = "; ".join("[" + ", ".join(p) + "]" for p in parts)
        raise ClaireError(
            f"Fragment '{fname}': the atoms are not all bonded to each other in {proj.ligand_pdb().name}; they form "
            f"{len(parts)} separate pieces: {desc}. Atom names come from your ligand model (molfile_to_params usually names "
            "atoms in file order, not by chemical numbering). Run `claire fragments --list-atoms` to see each atom's "
            "neighbours and choose a connected set.")


def _fragment_files(proj: Project, outdir: Path) -> list[tuple[str, Path]]:
    frags = proj.fragments()
    if not frags:
        raise ClaireError("No fragments defined. Add them under `fragments:` in claire.yaml "
                          "(atom names from `claire fragments --list-atoms`).")
    lig_lines = [l for l in proj.ligand_pdb().read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    obabel = None
    try:
        obabel = tools.executable(proj, "obabel")
    except ClaireError:
        print("[fragments] obabel not found; using RDKit to write mol files (bond orders from the ligand model).")
    out = []
    for i, f in enumerate(frags, 1):
        atoms = [a.strip() for a in f["atoms"]]
        sel = [l for l in lig_lines if l[12:16].strip() in atoms]
        missing = set(atoms) - {l[12:16].strip() for l in sel}
        if missing:
            raise ClaireError(f"Fragment {f.get('name', i)}: atoms {sorted(missing)} not in {proj.ligand_pdb().name}")
        _check_connected(proj, f.get("name", f"Fragment_{i}"), atoms)
        name = f"Fragment_{i}"
        pdb = outdir / f"{name}.pdb"
        pdb.write_text("\n".join(sel) + "\nEND\n")
        mol = outdir / f"{name}.mol"
        if obabel:
            subprocess.run([obabel, "-i", "pdb", pdb.name, "-o", "mol", "-O", mol.name], cwd=outdir,
                           capture_output=True, check=True)
        else:
            from rdkit import Chem
            full = Chem.MolFromMolFile(str(proj.stage("target", False) / f"{proj.lig}.sdf"), removeHs=False)
            if full is None:
                raise ClaireError("Install Open Babel (obabel) to convert fragments to mol files.")
            names = [l[12:16].strip() for l in lig_lines]
            keep = [names.index(a) for a in atoms]
            em = Chem.RWMol(full)
            for idx in sorted(set(range(full.GetNumAtoms())) - set(keep), reverse=True):
                em.RemoveAtom(idx)
            Chem.MolToMolFile(em.GetMol(), str(mol))
        out.append((name, mol))
    return out


def _rcsb_smiles(comp_ids: list[str]) -> dict:
    """Canonical CCD SMILES for chemical components, via the RCSB GraphQL API (batched)."""
    res = {}
    for i in range(0, len(comp_ids), 400):
        chunk = comp_ids[i:i + 400]
        q = '{chem_comps(comp_ids:%s){chem_comp{id} rcsb_chem_comp_descriptor{SMILES}}}' % json.dumps(chunk)
        req = urllib.request.Request("https://data.rcsb.org/graphql", data=json.dumps({"query": q}).encode(),
                                     headers={"Content-Type": "application/json", **sources.UA})
        data = json.loads(urllib.request.urlopen(req, timeout=120).read())["data"]["chem_comps"]
        for c in data:
            if c and c.get("rcsb_chem_comp_descriptor") and c["rcsb_chem_comp_descriptor"].get("SMILES"):
                res[c["chem_comp"]["id"]] = c["rcsb_chem_comp_descriptor"]["SMILES"]
    return res


def _rcsb_entries(comp_ids: list[str]) -> list[str]:
    entries = set()
    for i in range(0, len(comp_ids), 500):
        chunk = comp_ids[i:i + 500]
        entries.update(sources._search({
            "query": {"type": "terminal", "service": "text",
                      "parameters": {"attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                                     "operator": "in", "value": chunk}},
            "return_type": "entry", "request_options": {"return_all_hits": True}}))
    return sorted(entries)


def _rcsb_substructure(name: str, m) -> list[str]:
    """RCSB substructure search with a stereo-free SMILES (the relaxed graph match ignores stereo anyway),
    falling back to a Kekule SMILES if the service rejects the first form."""
    from rdkit import Chem
    import urllib.error
    forms = [Chem.MolToSmiles(m, isomericSmiles=False)]
    k = Chem.Mol(m)
    try:
        Chem.Kekulize(k, clearAromaticFlags=True)
        forms.append(Chem.MolToSmiles(k, isomericSmiles=False, kekuleSmiles=True))
    except Exception:
        pass
    last = None
    for smi in dict.fromkeys(forms):
        print(f"[fragments] {name}: searching RCSB for substructure {smi}")
        try:
            return sources.rcsb_ligands_with_substructure(smi)
        except urllib.error.HTTPError as e:
            last = e
        except urllib.error.URLError as e:
            raise ClaireError(f"Could not reach the RCSB search service ({e.reason}). Check your internet connection.")
    raise ClaireError(
        f"The RCSB search service returned HTTP {last.code} for {name} (SMILES {forms[0]}). This is usually a temporary "
        "server problem: wait a few minutes and run `claire fragments` again. If it persists, check that the fragment is "
        "a sensible substructure, or set sources.ligand_search: bundled to use the original offline database.")


def run(proj: Project, args) -> None:
    if getattr(args, "list_atoms", False):
        list_atoms(proj)
        return
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    out = proj.stage("fragments")
    frags = _fragment_files(proj, out)
    mode = proj.config["sources"].get("ligand_search", "rcsb")
    results = {}
    if mode == "bundled":
        smiles_db, lig2pdb = sources.bundled_db(proj)
    for name, molf in frags:
        m = Chem.MolFromMolFile(str(molf))
        if m is None:
            raise ClaireError(f"Could not read {molf}")
        if mode == "bundled":
            cands = {k: v["SMILES_ccd"] for k, v in smiles_db.items()}
        else:
            ids = _rcsb_substructure(name, m)
            print(f"[fragments] {name}: {len(ids)} candidate PDB ligands from RCSB substructure search")
            cands = _rcsb_smiles(ids)
        ligs = []
        for lig, smi in cands.items():                 # same check as the original Fragment_Search.py
            pm = Chem.MolFromSmiles(smi) if smi else None
            if pm is not None and pm.HasSubstructMatch(m):
                ligs.append(lig)
        if mode == "bundled":
            pdb_ids = sorted({p for l in ligs for p in lig2pdb.get(l, [])})
        else:
            pdb_ids = _rcsb_entries(ligs) if ligs else []
        results[name] = {"Ligands": ligs, "PDBs": pdb_ids, "SMILES": {l: cands[l] for l in ligs}}
        print(f"[fragments] {name}: {len(ligs)} ligands, {len(pdb_ids)} PDB entries")
    (out / "PDB_search_results.json").write_text(json.dumps(results))
    proj.write_manifest("fragments", {"search_results": "PDB_search_results.json", "mode": mode,
                                      "fragments": {n: {"mol": f"{n}.mol", "n_pdbs": len(results[n]["PDBs"])}
                                                    for n, _ in frags},
                                      "n_outputs": len(frags)})
    print("Next: `claire contacts`.")
