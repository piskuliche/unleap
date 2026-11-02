"""Fixtures: small systems built on the fly with RDKit + antechamber + tleap (needs AmberTools)."""

import os
import shutil
import subprocess

import pytest

if not os.environ.get("AMBERHOME") or shutil.which("tleap") is None:
    pytest.skip("AmberTools not available (source it first, e.g. `source $AMBERHOME/amber.sh`)", allow_module_level=True)

ACETAMINOPHEN = "CC(=O)Nc1ccc(O)cc1"
CHLOROBENZOATE = "[O-]C(=O)c1ccc(Cl)cc1"


def sh(cmd, cwd):
    res = subprocess.run(cmd, cwd=cwd, shell=True, capture_output=True, text=True)
    assert res.returncode == 0, f"{cmd}\n{res.stdout}\n{res.stderr}"
    return res.stdout


def tleap(script, cwd):
    with open(os.path.join(cwd, "leap.in"), "w") as fh:
        fh.write(script + "\nquit\n")
    out = sh("tleap -f leap.in", cwd)
    assert "Errors = 0" in out, out


def parametrize(smiles, cwd, atom_type="gaff2", charge=0):
    """SMILES -> LIG.mol2 + LIG.frcmod (antechamber gas charges for speed)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToMolFile(mol, os.path.join(cwd, "LIG.sdf"))
    sh(f"antechamber -i LIG.sdf -fi mdl -o LIG.mol2 -fo mol2 -c gas -nc {charge} "
       f"-at {atom_type} -rn LIG -s 0 -pf y", cwd)
    sh(f"parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod -s {atom_type}", cwd)


def build(cwd, preamble="source leaprc.gaff2", frcmods=("LIG.frcmod",), solute="sys = copy LIG",
          name="sys"):
    """Solvate + neutralize a system containing LIG; returns (parm7, rst7)."""
    loads = "\n".join(f"loadamberparams {f}" for f in frcmods)
    tleap(f"""{preamble}
source leaprc.water.tip3p
{loads}
LIG = loadmol2 LIG.mol2
{solute}
solvatebox sys TIP3PBOX 8
addions sys Na+ 0
addions sys Cl- 0
saveamberparm sys {name}.parm7 {name}.rst7""", cwd)
    return os.path.join(cwd, f"{name}.parm7"), os.path.join(cwd, f"{name}.rst7")


@pytest.fixture(scope="session")
def gaff2_box(tmp_path_factory):
    d = str(tmp_path_factory.mktemp("gaff2_box"))
    parametrize(ACETAMINOPHEN, d)
    return d, *build(d)


@pytest.fixture(scope="session")
def charged_complex(tmp_path_factory):
    """Charged ligand next to a capped peptide (ff14SB + GAFF2 + TIP3P)."""
    d = str(tmp_path_factory.mktemp("complex"))
    parametrize(CHLOROBENZOATE, d, charge=-1)
    return d, *build(d, preamble="source leaprc.protein.ff14SB\nsource leaprc.gaff2",
                     solute="pep = sequence {ACE ALA LYS ALA NME}\n"
                            "translate LIG {12 0 0}\nsys = combine {pep LIG}")
