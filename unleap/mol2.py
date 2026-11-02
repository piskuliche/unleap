"""TRIPOS mol2 writer: Amber atom types, exact charges, RDKit-inferred bond orders."""

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Geometry import Point3D


def bond_orders(lig):
    """Return ({(i, j): order} with i < j, kekulized, RDKit mol); (None, None) if RDKit cannot assign them."""
    mol = Chem.RWMol()
    for atom in lig.atoms:
        a = Chem.Atom(atom.atomic_number)
        a.SetNoImplicit(True)
        mol.AddAtom(a)
    for b in lig.bonds:
        mol.AddBond(b.atom1.idx, b.atom2.idx, Chem.BondType.SINGLE)
    conf = Chem.Conformer(len(lig.atoms))
    for atom in lig.atoms:
        conf.SetAtomPosition(atom.idx, Point3D(atom.xx, atom.xy, atom.xz))
    mol.AddConformer(conf, assignId=True)
    charge = int(round(sum(a.charge for a in lig.atoms)))
    try:
        rdDetermineBonds.DetermineBondOrders(mol, charge=charge, embedChiral=False)
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        return None, None
    orders = {}
    for b in mol.GetBonds():
        i, j = sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        orders[(i, j)] = int(b.GetBondTypeAsDouble())
    return orders, mol


def write_sdf(mol, ligand, path, name):
    """SDF with bond orders and formal charges: the RDKit/viewer-friendly copy (RDKit's mol2
    reader only understands SYBYL types, not GAFF)."""
    mol = Chem.Mol(mol)
    mol.SetProp("_Name", name)
    for atom, a in zip(mol.GetAtoms(), ligand.lig.atoms):
        atom.SetProp("molFileAlias", a.name)
    Chem.MolToMolFile(mol, path, kekulize=True)


def write(ligand, path, name, sdf_path=None):
    """Write the mol2 (and SDF, if bond orders were assigned); returns False when bond orders
    fell back to all-single."""
    lig = ligand.lig
    orders, mol = bond_orders(lig)
    if mol is not None and sdf_path:
        write_sdf(mol, ligand, sdf_path, name)
    bonds = sorted(tuple(sorted((b.atom1.idx, b.atom2.idx))) for b in lig.bonds)
    lines = ["@<TRIPOS>MOLECULE", name,
             f"{len(lig.atoms):5d} {len(bonds):5d} {1:5d} {0:5d} {0:5d}",
             "SMALL", "USER_CHARGES", "", "", "@<TRIPOS>ATOM"]
    for atom in lig.atoms:
        lines.append(f"{atom.idx + 1:7d} {atom.name:<4s} {atom.xx:14.4f}{atom.xy:11.4f}{atom.xz:11.4f} "
                     f"{atom.type:<8s} {1:4d} {name:<6s} {atom.charge:12.6f}")
    lines.append("@<TRIPOS>BOND")
    for n, (i, j) in enumerate(bonds, 1):
        order = orders.get((i, j), 1) if orders else 1
        lines.append(f"{n:6d}{i + 1:6d}{j + 1:6d} {order}")
    lines += ["@<TRIPOS>SUBSTRUCTURE",
              f"{1:6d} {name:<8s} {1:6d} TEMP              0 ****  ****    0 ROOT", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return orders is not None
