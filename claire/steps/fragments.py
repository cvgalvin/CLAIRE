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


def _fragment_files(proj: Project, outdir: Path) -> list[tuple[str, Path, str, str]]:
    """Write Fragment_<n>.pdb/.mol/.smarts/.search.smarts; returns (name, mol, search SMARTS, RCSB SMILES)."""
    from rdkit import Chem
    from .autofrag import fragment_queries, parent_mol
    frags = proj.fragments()
    if not frags:
        raise ClaireError("No fragments defined. Add them under `fragments:` in claire.yaml "
                          "(atom names from `claire fragments --list-atoms`), or use `fragments: auto`.")
    lig_lines = [l for l in proj.ligand_pdb().read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    mol, idx = parent_mol(proj)
    out = []
    for i, f in enumerate(frags, 1):
        atoms = [a.strip() for a in f["atoms"]]
        missing = set(atoms) - set(idx)
        if missing:
            raise ClaireError(f"Fragment {f.get('name', i)}: atoms {sorted(missing)} not in {proj.ligand_pdb().name}")
        _check_connected(proj, f.get("name", f"Fragment_{i}"), atoms)
        name = f"Fragment_{i}"
        sel = {a for a in atoms}
        (outdir / f"{name}.pdb").write_text("\n".join(l for l in lig_lines if l[12:16].strip() in sel) + "\nEND\n")
        # .mol in fragment atom order (coordinates for alignment); bonds from the whole ligand
        sub = Chem.RWMol()
        new = {}
        conf = Chem.Conformer(len(atoms))
        for k, a in enumerate(atoms):
            at = Chem.Atom(mol.GetAtomWithIdx(idx[a]).GetAtomicNum())
            new[a] = sub.AddAtom(at)
            conf.SetAtomPosition(k, mol.GetConformer().GetAtomPosition(idx[a]))
        kek = Chem.RWMol(mol)
        Chem.Kekulize(kek, clearAromaticFlags=True)
        for x in range(len(atoms)):
            for y in range(x + 1, len(atoms)):
                b = kek.GetBondBetweenAtoms(idx[atoms[x]], idx[atoms[y]])
                if b is not None:
                    sub.AddBond(new[atoms[x]], new[atoms[y]], b.GetBondType())
        sm = sub.GetMol()
        sm.AddConformer(conf, assignId=True)
        Chem.MolToMolFile(sm, str(outdir / f"{name}.mol"), kekulize=False)
        smarts_all, smarts_heavy, smiles = fragment_queries(proj, atoms)
        (outdir / f"{name}.smarts").write_text(smarts_all + "\n")          # extraction (incl. any listed H caps)
        (outdir / f"{name}.search.smarts").write_text(smarts_heavy + "\n")
        out.append((name, outdir / f"{name}.mol", smarts_heavy, smiles))
    return out


_CCD_MOLS = {}


def _ccd_mols(ccd: dict) -> dict:
    """RDKit molecules for all CCD SMILES (parsed once per process)."""
    if not _CCD_MOLS:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        for k, v in ccd.items():
            m = Chem.MolFromSmiles(v)
            if m is not None:
                _CCD_MOLS[k] = m
    return _CCD_MOLS


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


def _rcsb_substructure(name: str, m, smiles: str | None = None) -> list[str]:
    """RCSB substructure search with a stereo-free SMILES (the relaxed graph match ignores stereo anyway),
    falling back to a Kekule SMILES if the service rejects the first form."""
    from rdkit import Chem
    import urllib.error
    forms = [smiles] if smiles else [Chem.MolToSmiles(m, isomericSmiles=False)]
    try:
        k = Chem.Mol(m)
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


def show_auto(proj: Project, count: bool) -> None:
    import yaml
    from . import autofrag
    f = proj.stage("fragments") / "auto_fragments.yaml"
    if f.exists():
        f.unlink()                                  # recompute from the current params file
    d = proj.auto_fragments()
    print(f"Automatic fragments for {proj.lig} (written to {f.relative_to(proj.root)}):")
    for fr in d["fragments"]:
        kind = "H-bond group" if fr["name"] in d["hbond_fragments"] else "aromatic system"
        line = f"  {fr['name']:<22} {kind:<16} {', '.join(fr['atoms'])}"
        if count:
            from rdkit import Chem
            from .autofrag import fragment_queries
            _, smarts, smiles = fragment_queries(proj, fr["atoms"])
            q = Chem.MolFromSmarts(smarts)
            ok = [k for k, m in _ccd_mols(sources.ccd_smiles()).items() if m.HasSubstructMatch(q)]
            line += f"   [{len(ok)} PDB ligands]   {smarts}"
        print(line)
    print("\nSuggested polar_atoms (copy into claire.yaml and edit as needed):")
    print("\n".join("  " + l for l in yaml.safe_dump({"polar_atoms": d["polar_atoms_suggested"]},
                                                     default_flow_style=None, sort_keys=False).splitlines()))
    print("\nTo use these fragments set `fragments: auto` and `hbond_fragments: auto` in claire.yaml, then run "
          "`claire fragments`.")


def _fragment_mol(proj: Project, atoms):
    """RDKit mol for a set of ligand atoms, with bond orders from the params file."""
    from rdkit import Chem
    from .autofrag import chemistry
    _, heavy, el, nb, arom_bond, double = chemistry(proj)
    sel = [a for a in atoms if el.get(a) != "H"]
    m = Chem.RWMol()
    idx = {}
    for a in sel:
        at = Chem.Atom(el[a].capitalize())
        idx[a] = m.AddAtom(at)
    done = set()
    for a in sel:
        for b in nb[a]:
            if b in idx and frozenset((a, b)) not in done:
                done.add(frozenset((a, b)))
                bt = Chem.BondType.AROMATIC if arom_bond(a, b) else (Chem.BondType.DOUBLE if double(a, b) else Chem.BondType.SINGLE)
                m.AddBond(idx[a], idx[b], bt)
                if bt == Chem.BondType.AROMATIC:
                    m.GetAtomWithIdx(idx[a]).SetIsAromatic(True)
                    m.GetAtomWithIdx(idx[b]).SetIsAromatic(True)
    try:
        mol = m.GetMol()
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        return mol
    except Exception:
        return None


def run(proj: Project, args) -> None:
    if getattr(args, "list_atoms", False):
        list_atoms(proj)
        return
    if getattr(args, "auto", False):
        show_auto(proj, getattr(args, "count", False))
        return
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    out = proj.stage("fragments")
    frags = _fragment_files(proj, out)
    mode = proj.config["sources"].get("ligand_search", "ccd")
    results = {}
    if mode == "bundled":
        smiles_db, lig2pdb = sources.bundled_db(proj)
    if mode == "ccd":
        ccd = sources.ccd_smiles()
        ccd_mols = _ccd_mols(ccd)
        print(f"[fragments] searching {len(ccd_mols)} chemical components of the current CCD locally")
    for name, molf, smarts, smiles in frags:
        q = Chem.MolFromSmarts(smarts)
        if mode == "ccd":
            ligs = [k for k, m in ccd_mols.items() if m.HasSubstructMatch(q)]
            pdb_ids = _rcsb_entries(ligs) if ligs else []
            results[name] = {"Ligands": ligs, "PDBs": pdb_ids, "SMILES": {l: ccd[l] for l in ligs}}
            print(f"[fragments] {name}: {len(ligs)} ligands, {len(pdb_ids)} PDB entries")
            continue
        if mode == "bundled":
            cands = {k: v["SMILES_ccd"] for k, v in smiles_db.items()}
        else:
            ids = []
            for form in smiles.split("|"):        # Kekule and aromatic forms; union of candidates
                ids += [x for x in _rcsb_substructure(name, Chem.MolFromSmiles(form, sanitize=False) or Chem.MolFromSmarts(smarts), form)
                        if x not in ids]
            print(f"[fragments] {name}: {len(ids)} candidate PDB ligands from RCSB substructure search")
            cands = _rcsb_smiles(ids)
        ligs = []
        for lig, smi in cands.items():          # charge-agnostic match; aromaticity from the whole target ligand
            pm = Chem.MolFromSmiles(smi) if smi else None
            if pm is not None and pm.HasSubstructMatch(q):
                ligs.append(lig)
        if mode == "bundled":
            pdb_ids = sorted({p for l in ligs for p in lig2pdb.get(l, [])})
        else:
            pdb_ids = _rcsb_entries(ligs) if ligs else []
        results[name] = {"Ligands": ligs, "PDBs": pdb_ids, "SMILES": {l: cands[l] for l in ligs}}
        print(f"[fragments] {name}: {len(ligs)} ligands, {len(pdb_ids)} PDB entries")
    (out / "PDB_search_results.json").write_text(json.dumps(results))
    proj.write_manifest("fragments", {"search_results": "PDB_search_results.json", "mode": mode,
                                      "fragments": {n: {"mol": f"{n}.mol", "smarts": f"{n}.smarts", "query_smiles": smi,
                                                        "n_ligands": len(results[n]["Ligands"]), "n_pdbs": len(results[n]["PDBs"])}
                                                    for n, _, _, smi in frags},
                                      "n_outputs": len(frags)})
    print("Next: `claire contacts`.")
