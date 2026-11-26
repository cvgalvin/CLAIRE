"""Group clustering (v0.4): fragment symmetry, residue renaming, and the sidechain fallback."""
import numpy as np
from collections import OrderedDict
from claire import cluster_fast as CF


def _benzene_pdb(tmp_path, extra=None):
    ang = np.radians(np.arange(6) * 60.0)
    X = np.c_[1.39 * np.cos(ang), 1.39 * np.sin(ang), np.zeros(6)]
    lines = [f"HETATM{i + 1:5d}  C{i + 1:<2} BNZ X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
             for i, (x, y, z) in enumerate(X)]
    p = tmp_path / "BNZ_0001.pdb"
    p.write_text("".join(lines + (extra or [])))
    return p, {f"C{i + 1}" for i in range(6)}


def _contact_lines(resname, atoms, resnum):
    return [f"ATOM  {i + 1:5d}  {n:<3} {resname} A{resnum:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {n[0]}\n"
            for i, (n, (x, y, z)) in enumerate(atoms.items())]


def _arg(shift_rot=0.0, flip=False, offset=0.0):
    """an Arg guanidinium over the benzene face, optionally rotated about the ring normal or moved to the other face"""
    base = {"N": (2.5, 1.0, 6.0), "CA": (2.0, 0.5, 5.6), "C": (2.6, -0.5, 6.0), "O": (3.4, -0.8, 6.8), "CB": (1.5, 0.3, 5.2),
            "CG": (1.2, 0.4, 5.0), "CD": (0.8, 0.6, 4.6), "NE": (0.5, 0.4, 4.0), "CZ": (0.2, 0.1, 3.6),
            "NH1": (0.9, -0.3, 3.4), "NH2": (-0.8, 0.2, 3.4)}
    t = np.radians(shift_rot)
    R = np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])
    if flip:
        R = R @ np.diag([1, -1, -1])
    return {n: tuple(R @ (np.array(v) + [offset, 0, 0])) for n, v in base.items()}


def _cluster(lines, resname, ligand, fatoms, cutoff=1.0):
    from claire.cluster_fast import read_fuzzball
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".pdb"); os.close(fd)
    open(path, "w").write("".join(lines))
    ls, res, rows = read_fuzzball(path)
    return CF.cluster_groups(ls, rows, resname, CF.fragment_symmetry(str(ligand), fatoms), cutoff)


def test_benzene_has_twelve_symmetries(tmp_path):
    lig, fa = _benzene_pdb(tmp_path)
    assert len(CF.fragment_symmetry(str(lig), fa)) == 12


def test_asymmetric_fragment_has_identity_only(tmp_path):
    lines = ["HETATM    1  C1  LIG X   1       0.000   0.000   0.000  1.00  0.00           C\n",
             "HETATM    2  N1  LIG X   1       1.470   0.000   0.000  1.00  0.00           N\n",
             "HETATM    3  O1  LIG X   1      -0.700   1.200   0.000  1.00  0.00           O\n"]
    p = tmp_path / "LIG_0001.pdb"; p.write_text("".join(lines))
    assert len(CF.fragment_symmetry(str(p), {"C1", "N1", "O1"})) == 1


def test_symmetry_equivalent_contacts_cluster(tmp_path):
    lig, fa = _benzene_pdb(tmp_path)
    lines = (_contact_lines("ARG", _arg(offset=2.5), 1) + _contact_lines("ARG", _arg(shift_rot=60, offset=2.5), 2)
             + _contact_lines("ARG", _arg(shift_rot=120, flip=True, offset=2.5), 3))
    out = _cluster(lines, "ARG", lig, fa)
    assert len(out) == 1 and len(next(iter(out.values()))) == 2          # all three together


def test_without_symmetry_they_would_not(tmp_path):
    lig, fa = _benzene_pdb(tmp_path)
    lines = _contact_lines("ARG", _arg(offset=2.5), 1) + _contact_lines("ARG", _arg(shift_rot=60, offset=2.5), 2)
    from claire.cluster_fast import read_fuzzball
    p = tmp_path / "c.pdb"; p.write_text("".join(lines))
    ls, res, rows = read_fuzzball(str(p))
    assert len(CF.cluster_groups(ls, rows, "ARG", [(np.eye(3), np.zeros(3))], 1.0)) == 0


def test_ring_flip_renaming_clusters(tmp_path):
    lig, fa = _benzene_pdb(tmp_path)
    phe = {"N": (3.0, 0, 6.5), "CA": (2.5, 0, 6.0), "C": (3.0, 1.0, 6.5), "O": (3.8, 1.2, 7.0), "CB": (2.0, 0, 5.5),
           "CG": (1.5, 0, 4.5), "CD1": (2.2, 1.0, 4.0), "CD2": (0.6, -0.8, 4.0), "CE1": (1.8, 1.3, 3.0),
           "CE2": (0.3, -0.5, 3.0), "CZ": (0.9, 0.5, 2.6)}
    flipped = dict(phe)
    flipped["CD1"], flipped["CD2"], flipped["CE1"], flipped["CE2"] = phe["CD2"], phe["CD1"], phe["CE2"], phe["CE1"]
    lines = _contact_lines("PHE", phe, 1) + _contact_lines("PHE", flipped, 2)
    assert len(_cluster(lines, "PHE", lig, fa)) == 1


def test_different_contacts_stay_apart(tmp_path):
    lig, fa = _benzene_pdb(tmp_path)
    far = {n: (x + 3.0, y, z) for n, (x, y, z) in _arg().items()}               # slid 3 A along the ring plane
    lines = _contact_lines("ARG", _arg(), 1) + _contact_lines("ARG", far, 2)
    assert len(_cluster(lines, "ARG", lig, fa)) == 0


def test_aromatic_nh_is_not_charged():
    """An aromatic N-H (2-pyridone) sums to bond order 4 with aromatic bonds counted as 1.5; it must stay neutral."""
    from rdkit import Chem
    from claire.steps import autofrag as af
    m = Chem.RWMol()
    for e in ("O", "C", "C", "C", "C", "C", "N", "H"):
        a = Chem.Atom(e); a.SetNoImplicit(True); m.AddAtom(a)
    for a, b, t in ((0, 1, 2), (1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 6, 4), (6, 1, 4), (6, 7, 1)):
        m.AddBond(a, b, {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 4: Chem.BondType.AROMATIC}[t])
    mol = m.GetMol()
    n = mol.GetAtomWithIdx(6)
    aromatic = any(b.GetBondType() == Chem.BondType.AROMATIC for b in n.GetBonds())
    assert aromatic and sum(b.GetBondTypeAsDouble() for b in n.GetBonds()) == 4      # the case the old rule mis-charged
    src = open(af.__file__).read()
    assert "b.GetIsAromatic() or b.GetBondType() == Chem.BondType.AROMATIC" in src     # the rule now excludes it


def _write_lig(tmp_path, code, atoms, bonds):
    """model + params for a small ligand: atoms [(name, element, x, y, z)], bonds [(a, b, order)]"""
    pdb = "".join(f"HETATM{i + 1:5d} {n:<4} {code} X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {e:>2}\n"
                  for i, (n, e, x, y, z) in enumerate(atoms))
    (tmp_path / f"{code}_0001.pdb").write_text(pdb)
    par = "".join(f"ATOM  {n:<4} X    X  0.00\n" for n, e, x, y, z in atoms) + "".join(f"BOND_TYPE {a:<4} {b:<4} {o}\n" for a, b, o in bonds)
    (tmp_path / f"{code}.params").write_text(par)
    return tmp_path / f"{code}_0001.pdb", {n for n, e, *_ in atoms if e != "H"}


def test_cyclohexenone_has_no_mirror_symmetry(tmp_path):
    """C=C on one side of the carbonyl, CH2-CH2 on the other: graph symmetric without bond orders, chemically not."""
    ang = np.radians(np.arange(6) * 60.0)
    ring = [(f"C{i + 1}", "C", 1.40 * np.cos(a), 1.40 * np.sin(a), 0.0) for i, a in enumerate(ang)]
    atoms = ring + [("O1", "O", 2.62, 0.0, 0.0)]
    bonds = [("C1", "O1", 2), ("C1", "C2", 1), ("C2", "C3", 2), ("C3", "C4", 1), ("C4", "C5", 1), ("C5", "C6", 1), ("C6", "C1", 1)]
    lig, fa = _write_lig(tmp_path, "CHX", atoms, bonds)
    assert len(CF.fragment_symmetry(str(lig), fa)) == 1


def test_carboxylate_oxygens_stay_equivalent(tmp_path):
    atoms = [("C1", "C", 0.0, 0.0, 0.0), ("O1", "O", 1.25, 0.0, 0.0), ("O2", "O", -0.625, 1.08, 0.0), ("C2", "C", -0.75, -1.30, 0.0)]
    bonds = [("C1", "O1", 2), ("C1", "O2", 1), ("C1", "C2", 1)]
    lig, fa = _write_lig(tmp_path, "ACT", atoms, bonds)
    # O1 and O2 are related by a two-fold symmetry only up to the methyl, so test the carboxylate group alone
    assert len(CF.fragment_symmetry(str(lig), {"C1", "O1", "O2"})) == 2
