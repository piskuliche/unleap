# unleap

[![tests](https://github.com/piskuliche/unleap/actions/workflows/tests.yml/badge.svg)](https://github.com/piskuliche/unleap/actions/workflows/tests.yml)
[![docs](https://github.com/piskuliche/unleap/actions/workflows/docs.yml/badge.svg)](https://piskuliche.github.io/unleap/)
[![PyPI](https://img.shields.io/pypi/v/unleap.svg)](https://pypi.org/project/unleap/)

**Get a ligand's parameter files back out of a built Amber system.**

tleap turns a ligand's `mol2` / `frcmod` / `lib` into a system (`parm7` + `rst7`); `unleap` goes
the other way. It strips the ligand out of a complex or solvated system and regenerates its
`frcmod` (a diff against GAFF/GAFF2, plus a self-contained copy), OFF `lib`, `mol2`, and an `sdf`
for RDKit. It then rebuilds the ligand with tleap to prove the files reproduce the original
parameters and energies.

## Install

unleap needs AmberTools (tleap, sander), which is not pip-installable. The simplest setup:

```bash
conda create -n unleap -c conda-forge ambertools rdkit
conda activate unleap
pip install unleap
```

## Use

```bash
unleap -p complex.parm7 -c complex.rst7 -o lig_out                      # ligand auto-detected
unleap -p unisc.parm7 -c 0.00000000_ti.rst7 --resname L00 -o L00_out   # RBFE dual topology
unleap -p complex.parm7 -c traj.nc --frame -1 -o lig_out               # from a trajectory
```

Hydrogen mass repartitioning is detected and undone. The output directory gets `LIG.frcmod`,
`LIG.full.frcmod`, `LIG.lib`, `LIG.mol2`, `LIG.sdf`, a `tleap.in` that loads them, and a
`report.txt` with the PASS/FAIL validation results.

## Documentation

**<https://piskuliche.github.io/unleap/>**: the user guide, developer guide and API reference.

To build it locally: `pip install -e ".[docs]" && make -C docs html`.

## Development

```bash
git clone https://github.com/piskuliche/unleap.git && cd unleap
pip install -e ".[test,docs]"
python -m pytest -q          # needs AmberTools on PATH; skipped otherwise
```

## License

MIT
