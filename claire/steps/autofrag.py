"""Automatic fragment definition from the ligand's chemistry (atom names and bond orders from the params file).

Fragments:
  * one per hydrogen-bonding group: the polar atom, all heavy atoms within two bonds, double-bond partners of those
    atoms, and any aromatic ring system they touch (so a phenol gives ring + OH, a ketone gives the carbonyl with its
    substituents and any conjugated C=C). Groups with identical atom sets (e.g. both carboxylate oxygens) are merged.
  * one per (fused) aromatic ring system.
Suggested polar_atoms entries are produced for every acceptor and donor (donors list the polar hydrogen first).
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from ..project import ClaireError, Project


def read_params(params: Path):
    atoms, bonds = [], []
    for line in Path(params).read_text().splitlines():
        if line.startswith("ATOM "):
            atoms.append(line.split()[1])
        elif line.startswith("BOND"):
            parts = line.split()
            order = parts[3] if len(parts) > 3 and line.startswith("BOND_TYPE") else "1"
            bonds.append((parts[1], parts[2], order))
    return atoms, bonds


def _elements(proj: Project) -> dict:
    el = {}
    for l in proj.ligand_pdb().read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            n = l[12:16].strip()
            el[n] = (l[76:78].strip() or "".join(c for c in n if c.isalpha())[:1]).upper()
    return el


def chemistry(proj: Project):
    atoms, bonds = read_params(proj.params())
    el = _elements(proj)
    missing = [a for a in atoms if a not in el]
    if missing:
        raise ClaireError(f"Atoms {missing[:5]} from the params file are not in {proj.ligand_pdb().name}.")
    nb = defaultdict(list)
    order = {}
    for a, b, o in bonds:
        nb[a].append(b)
        nb[b].append(a)
        order[frozenset((a, b))] = o
    heavy = [a for a in atoms if el[a] != "H"]
    perceived = _perceive_orders(proj, atoms, bonds, el)
    if perceived is not None:
        order = perceived
        mismatch = [tuple(sorted(k)) for k, v in perceived.items()
                    if v != "4" and {"1": "1", "2": "2", "3": "3"}.get(next((o for a, b, o in bonds if frozenset((a, b)) == k), "1"), "1") != v]
        if mismatch:
            print(f"[fragments] note: bond orders in {proj.params().name} differ from those implied by the model's "
                  f"geometry and hydrogens for {len(mismatch)} bond(s), e.g. {mismatch[:3]}; using the perceived orders.")
    arom_bond = lambda a, b: order.get(frozenset((a, b)), "1").upper() in ("4", "ARO", "AR", "A")
    double = lambda a, b: order.get(frozenset((a, b)), "1") in ("2", "3")
    return atoms, heavy, el, nb, arom_bond, double


def parent_mol(proj: Project):
    """Whole ligand as an RDKit molecule (atom order = params order, explicit H, perceived bond orders and aromaticity,
    model coordinates); atom names in the 'name' property. Returns (mol, {name: index})."""
    from rdkit import Chem
    atoms, bonds = read_params(proj.params())
    el = _elements(proj)
    orders = _perceive_orders(proj, atoms, bonds, el)
    m = Chem.RWMol()
    idx = {}
    for a in atoms:
        at = Chem.Atom(el[a].capitalize())
        at.SetNoImplicit(True)
        at.SetProp("name", a)
        idx[a] = m.AddAtom(at)
    bt = {"1": Chem.BondType.SINGLE, "2": Chem.BondType.DOUBLE, "3": Chem.BondType.TRIPLE, "4": Chem.BondType.AROMATIC}
    for a, b, o in bonds:
        k = frozenset((a, b))
        m.AddBond(idx[a], idx[b], bt.get((orders or {}).get(k, o), Chem.BondType.SINGLE))
    xyz = {l[12:16].strip(): (float(l[30:38]), float(l[38:46]), float(l[46:54]))
           for l in proj.ligand_pdb().read_text().splitlines() if l.startswith(("ATOM", "HETATM"))}
    conf = Chem.Conformer(len(atoms))
    for a, i in idx.items():
        conf.SetAtomPosition(i, xyz[a])
    mol = m.GetMol()
    mol.AddConformer(conf, assignId=True)
    for a in mol.GetAtoms():                              # formal charges implied by valence (e.g. sulfonate O-)
        if a.GetSymbol() == "O" and a.GetDegree() == 1 and mol.GetBondBetweenAtoms(a.GetIdx(), a.GetNeighbors()[0].GetIdx()).GetBondType() == Chem.BondType.SINGLE \
                and a.GetNeighbors()[0].GetSymbol() != "H":
            a.SetFormalCharge(-1)
        if a.GetSymbol() == "N" and sum(b.GetBondTypeAsDouble() for b in a.GetBonds()) == 4:
            a.SetFormalCharge(1)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    return mol, idx


def fragment_queries(proj: Project, atoms: list[str]):
    """(extraction SMARTS incl. listed H, search SMARTS heavy atoms only, RCSB query SMILES) for a fragment.
    Aromaticity comes from the whole ligand; charges are left unspecified; atom map numbers give file order."""
    from rdkit import Chem
    mol, idx = parent_mol(proj)
    q = Chem.RWMol(mol)
    for a in q.GetAtoms():
        a.SetFormalCharge(0)
        a.SetNoImplicit(False)
    order = [idx[a] for a in atoms]
    heavy = [i for i in order if q.GetAtomWithIdx(i).GetAtomicNum() != 1]

    import re
    # polar heavy atoms that carry hydrogens in the ligand keep their H count in the query (OH stays OH, so a phenol
    # fragment does not match aryl ethers or esters)
    hcount = {i: sum(1 for nb_ in mol.GetAtomWithIdx(i).GetNeighbors() if nb_.GetAtomicNum() == 1)
              for i in order if mol.GetAtomWithIdx(i).GetSymbol() in ("O", "N", "S")}

    def to_smarts(use, maps):
        for a in q.GetAtoms():
            a.SetAtomMapNum(0)
        for pos, i in enumerate(use, 1):
            q.GetAtomWithIdx(i).SetAtomMapNum(pos)
        smi = Chem.MolFragmentToSmiles(q, atomsToUse=use, isomericSmiles=False, canonical=False, kekuleSmiles=False)

        def tok(m):
            body, n = m.group(1), int(m.group(2))
            i = use[n - 1]
            if body == "H":
                body = "#1"
            h = hcount.get(i, 0)
            if h and "H" not in body.lstrip("#"):
                body = f"{body};H{h}"
            return f"[{body}:{n}]" if maps else (f"[{body}]" if (";" in body or len(body) > 2 or body.startswith("#")) else body)
        return re.sub(r"\[([^\]:]+):(\d+)\]", tok, smi)

    smarts_all = to_smarts(order, True)
    smarts_heavy = to_smarts(heavy, False)
    k = Chem.RWMol(mol)
    for a in k.GetAtoms():
        a.SetFormalCharge(0)
    Chem.Kekulize(k, clearAromaticFlags=True)
    smiles = Chem.MolFragmentToSmiles(k, atomsToUse=heavy, kekuleSmiles=True, isomericSmiles=False)
    # RCSB's substructure search is sensitive to the bond representation: also offer the aromatic form when valid
    arom = Chem.MolFragmentToSmiles(q, atomsToUse=heavy, isomericSmiles=False)
    if Chem.MolFromSmiles(arom) is not None and arom != smiles:
        smiles = smiles + "|" + arom
    return smarts_all, smarts_heavy, smiles


def _perceive_orders(proj: Project, atoms, bonds, el):
    """Bond orders from connectivity + explicit hydrogens (RDKit), robust to wrong orders in the params file."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError:
        return None
    m = Chem.RWMol()
    idx = {}
    for a in atoms:
        at = Chem.Atom(el[a].capitalize())
        at.SetNoImplicit(True)
        idx[a] = m.AddAtom(at)
    for a, b, _ in bonds:
        m.AddBond(idx[a], idx[b], Chem.BondType.SINGLE)
    xyz = {l[12:16].strip(): (float(l[30:38]), float(l[38:46]), float(l[46:54]))
           for l in proj.ligand_pdb().read_text().splitlines() if l.startswith(("ATOM", "HETATM"))}
    conf = Chem.Conformer(len(atoms))
    for a, i in idx.items():
        conf.SetAtomPosition(i, xyz[a])
    m.AddConformer(conf, assignId=True)
    for charge in (0, -1, 1, -2, 2):
        mol = Chem.Mol(m)
        try:
            rdDetermineBonds.DetermineBondOrders(mol, charge=charge)
            Chem.SanitizeMol(mol)
        except Exception:
            continue
        name = {i: a for a, i in idx.items()}
        out = {}
        for bd in mol.GetBonds():
            k = frozenset((name[bd.GetBeginAtomIdx()], name[bd.GetEndAtomIdx()]))
            out[k] = "4" if bd.GetIsAromatic() else {Chem.BondType.SINGLE: "1", Chem.BondType.DOUBLE: "2",
                                                      Chem.BondType.TRIPLE: "3"}.get(bd.GetBondType(), "1")
        return out
    return None


def _aromatic_systems(heavy, nb, arom_bond):
    """Connected components of aromatic bonds (fused ring systems)."""
    seen, systems = set(), []
    for a in heavy:
        if a in seen or not any(arom_bond(a, b) for b in nb[a]):
            continue
        comp, stack = set(), [a]
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack += [y for y in nb[x] if arom_bond(x, y) and y not in comp]
        seen |= comp
        if len(comp) >= 5:
            systems.append(comp)
    return systems


def _small_rings(heavy, hnb, max_size=6):
    """Simple cycles of up to max_size atoms (as frozensets)."""
    rings = set()
    for start in heavy:
        stack = [(start, [start])]
        while stack:
            x, path = stack.pop()
            for y in hnb(x):
                if y == start and len(path) >= 3:
                    rings.add(frozenset(path))
                elif y not in path and len(path) < max_size:
                    stack.append((y, path + [y]))
    return [set(r) for r in rings]


def polar_roles(heavy, el, nb, arom_bond, double):
    """{atom: 'acc'|'don'|'accdon'} for ligand H-bond donors/acceptors."""
    roles = {}
    for a in heavy:
        e = el[a]
        hs = [x for x in nb[a] if el[x] == "H"]
        hv = [x for x in nb[a] if el[x] != "H"]
        if e == "O":
            roles[a] = "accdon" if hs else "acc"
        elif e == "N":
            conj = any(double(a, x) or arom_bond(a, x) or any(double(x, y) or arom_bond(x, y) for y in nb[x] if y != a)
                       for x in hv)
            aromatic = any(arom_bond(a, x) for x in hv)
            acc = (aromatic and len(hv) == 2 and not hs) or (not conj and not aromatic and len(hv) + len(hs) <= 3) \
                or any(double(a, x) for x in hv)
            don = bool(hs)
            if acc or don:
                roles[a] = ("acc" if acc else "") + ("don" if don else "")
    return roles


def define(proj: Project) -> dict:
    atoms, heavy, el, nb, arom_bond, double = chemistry(proj)
    systems = _aromatic_systems(heavy, nb, arom_bond)
    roles = polar_roles(heavy, el, nb, arom_bond, double)
    hnb = lambda a: [x for x in nb[a] if el[x] != "H"]
    from collections import deque

    def bond_dist(a, b):
        seen, dq = {a: 0}, deque([a])
        while dq:
            x = dq.popleft()
            if x == b:
                return seen[x]
            for y in hnb(x):
                if y not in seen:
                    seen[y] = seen[x] + 1
                    dq.append(y)
        return 99

    rings = _small_rings(heavy, hnb)

    def group_for(a):
        """Polar atom + atoms within two bonds; multiple-bond partners; aromatic systems bonded directly to the polar
        atom or its neighbours' multiple-bond partners; small non-aromatic rings sharing >= 2 atoms with that core;
        exocyclic double-bond partners. No cascading through fused rings."""
        first = {a} | set(hnb(a))
        core = set(first)
        for x in hnb(a):
            core.update(hnb(x))
        g = set(core)
        for x in list(core):
            g.update(y for y in hnb(x) if double(x, y))
        shell = first | {y for x in first for y in hnb(x) if double(x, y)}
        for s_ in systems:
            if s_ & shell:
                g |= s_
        added = []
        for r in rings:
            if len(r & core) >= 2 and not any(r <= s_ for s_ in systems):
                g |= r
                added.append(r)
        for r in added:                       # aromatic systems fused to an added ring (one level, no cascade)
            for s_ in systems:
                if len(r & s_) >= 2:
                    g |= s_
        for x in list(g):
            g.update(y for y in hnb(x) if double(x, y))
        return g

    groups = {}
    for a in roles:
        groups.setdefault(frozenset(group_for(a)), []).append(a)
    # merge groups whose polar atoms are within three bonds (vicinal diols, catechols, carboxylates, lactones)
    merged = True
    while merged:
        merged = False
        items = list(groups.items())
        for i, (g1, p1) in enumerate(items):
            for g2, p2 in items[i + 1:]:
                if min(bond_dist(x, y) for x in p1 for y in p2) <= 3:
                    groups.pop(g1)
                    groups.pop(g2)
                    groups[frozenset(g1 | g2)] = sorted(set(p1) | set(p2))
                    merged = True
                    break
            if merged:
                break
    frags, hb_names = [], []
    for g, polar in groups.items():
        name = "hb_" + "_".join(sorted(polar))
        frags.append({"name": name, "atoms": sorted(g, key=atoms.index), "polar": sorted(polar)})
        hb_names.append(name)
    for i, s in enumerate(systems, 1):
        # aromatic systems already contained in a larger (H-bond group) fragment are not run separately
        if not any(set(s) <= set(f["atoms"]) for f in frags):
            frags.append({"name": f"aromatic_{i}", "atoms": sorted(s, key=atoms.index)})
    # suggested polar_atoms: [polar atom, bonded heavy atom, next heavy atom]; donors start with the polar H
    suggestions = []
    for a, role in roles.items():
        b = next((x for x in hnb(a)), None)
        c = next((y for y in (hnb(b) if b else []) if y != a), None)
        if b is None or c is None:
            continue
        conj = any(double(a, x) or arom_bond(a, x) for x in hnb(a)) or any(arom_bond(b, y) for y in hnb(b))
        hyb = "sp2" if conj else "sp3"     # carbonyls, aromatic N, and O/N attached to an aromatic ring (e.g. phenol)
        if "acc" in role:
            suggestions.append({"atoms": [a, b, c], "donor": False, "hybridization": hyb})
        if "don" in role:
            h = next(x for x in nb[a] if el[x] == "H")
            suggestions.append({"atoms": [h, a, b], "donor": True, "hybridization": hyb})
    return {"fragments": frags, "hbond_fragments": hb_names, "polar_atoms_suggested": suggestions, "roles": roles}
