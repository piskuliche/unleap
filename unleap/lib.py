"""Write the ligand as an Amber OFF library unit."""

from parmed.amber import AmberOFFLibrary
from parmed.modeller import ResidueTemplate


def write(ligand, path, name):
    tmpl = ResidueTemplate.from_residue(ligand.lig.residues[0])
    tmpl.name = name
    tmpl.head = tmpl.tail = None
    AmberOFFLibrary.write({name: tmpl}, path)
