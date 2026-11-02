"""Write the ligand as an Amber OFF library unit."""

import io

from parmed.amber import AmberOFFLibrary
from parmed.modeller import ResidueTemplate


def render(ligand, name):
    """OFF library text holding one unit, ``name``."""
    tmpl = ResidueTemplate.from_residue(ligand.lig.residues[0])
    tmpl.name = name
    tmpl.head = tmpl.tail = None
    buf = io.StringIO()
    AmberOFFLibrary.write({name: tmpl}, buf)
    return buf.getvalue()


def write(ligand, path, name):
    with open(path, "w") as fh:
        fh.write(render(ligand, name))
