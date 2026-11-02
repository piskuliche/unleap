import os
import re

import parmed as pmd
import pytest

from conftest import ACETAMINOPHEN, build, parametrize, tleap
from unleap.cli import main


def run(parm, rst, outdir, *extra):
    rc = main(["-p", parm, "-c", rst, "-o", str(outdir), *extra])
    report = open(os.path.join(outdir, "report.txt")).read() if os.path.exists(
        os.path.join(outdir, "report.txt")) else ""
    return rc, report


def sections(path):
    """Parse a frcmod into {section: [lines]}."""
    out, cur = {}, None
    for line in open(path).read().splitlines()[1:]:
        key = line.strip().upper()
        if key in ("MASS", "BOND", "ANGLE", "ANGL", "DIHE", "IMPROPER", "IMPR", "NONB", "NONBON"):
            cur = key[:4]
            out.setdefault(cur, [])
        elif line.strip() and cur:
            out[cur].append(line)
    return out


def dihe_keys(path, sec="DIHE"):
    keys = set()
    for line in sections(path).get(sec, []):
        quad = tuple(t.strip() for t in line[:11].split("-"))
        keys.add(min(quad, quad[::-1]) if sec == "DIHE" else quad)
    return keys


def test_gaff2_roundtrip(gaff2_box, tmp_path):
    d, parm, rst = gaff2_box
    rc, report = run(parm, rst, tmp_path)
    assert rc == 0, report
    assert "base force field: gaff2" in report
    assert report.count("[PASS]") == 3
    for ext in (".mol2", ".frcmod", ".full.frcmod", ".lib"):
        assert os.path.getsize(tmp_path / f"LIG{ext}") > 0
    # every dihedral the diff needs, parmchk2 also had to supply
    assert dihe_keys(tmp_path / "LIG.frcmod") <= dihe_keys(os.path.join(d, "LIG.frcmod"))
    # charges are exact
    orig = pmd.load_file(os.path.join(d, "LIG.mol2"))
    new = pmd.load_file(str(tmp_path / "LIG.mol2"))
    assert [round(a.charge, 6) for a in orig.atoms] == [round(a.charge, 6) for a in new.atoms]


def test_custom_parameters_are_recovered(gaff2_box, tmp_path_factory, tmp_path):
    src, _, _ = gaff2_box
    d = str(tmp_path_factory.mktemp("custom"))
    for f in ("LIG.mol2", "LIG.frcmod"):
        os.symlink(os.path.join(src, f), os.path.join(d, f))
    with open(os.path.join(d, "custom.frcmod"), "w") as fh:
        fh.write("custom edits\nMASS\noh 16.00\n\nDIHE\nc3-c -ns-ca   1    3.000       0.000           2.000\n\n"
                 "NONBON\n  oh          1.9000  0.2000\n\n")
    parm, rst = build(d, frcmods=("LIG.frcmod", "custom.frcmod"))
    rc, report = run(parm, rst, tmp_path)
    assert rc == 0, report
    sec = sections(tmp_path / "LIG.frcmod")
    dihe = [l for l in sec["DIHE"] if l.startswith("c3-c -ns-ca")]
    assert len(dihe) == 1 and "3.00000000" in dihe[0] and " 0.000 " in dihe[0]
    assert any(l.startswith("oh") and "1.90000000" in l and "0.20000000" in l for l in sec["NONB"])
    assert [l.split()[0] for l in sec["NONB"]] == ["oh"]


def test_complex_autodetect_and_bond_orders(charged_complex, tmp_path):
    d, parm, rst = charged_complex
    rc, report = run(parm, rst, tmp_path)
    assert rc == 0, report
    assert "residue LIG" in report and "net charge -1.0000" in report
    bonds = open(tmp_path / "LIG.mol2").read().split("@<TRIPOS>BOND")[1].split("@<TRIPOS>")[0]
    orders = [line.split()[3] for line in bonds.strip().splitlines()]
    assert orders.count("2") == 4      # 3 ring + 1 carboxylate C=O (kekulized)
    lib = open(tmp_path / "LIG.lib").read()
    assert "ALA" not in lib and "WAT" not in lib


def test_hmr_is_undone(gaff2_box, tmp_path_factory, tmp_path):
    from parmed.tools import HMassRepartition
    d, parm, rst = gaff2_box
    hmr = pmd.load_file(parm, rst)
    HMassRepartition(hmr).execute()
    hdir = tmp_path_factory.mktemp("hmr")
    hparm = str(hdir / "hmr.parm7")
    hmr.save(hparm)
    rc, report = run(hparm, rst, tmp_path)
    assert rc == 0, report
    assert "HMR detected" in report
    masses = {l.split()[0]: float(l.split()[1]) for l in sections(tmp_path / "LIG.full.frcmod")["MASS"]}
    orig = pmd.load_file(parm)
    for a in orig.view[":LIG"].atoms:
        assert masses[a.type] == pytest.approx(a.mass, abs=1e-4), a.type


def test_keep_masses_rejects_split_types(gaff2_box, tmp_path_factory, tmp_path, capsys):
    from parmed.tools import HMassRepartition
    d, parm, rst = gaff2_box
    hmr = pmd.load_file(parm, rst)
    HMassRepartition(hmr).execute()
    hparm = str(tmp_path_factory.mktemp("hmr_keep") / "hmr.parm7")
    hmr.save(hparm)
    # HMR gives ca-with-H and ca-without-H different masses: not expressible as one frcmod type
    rc, _ = run(hparm, rst, tmp_path, "--keep-masses", "--no-check")
    assert rc == 2 and "hydrogen mass repartitioning" in capsys.readouterr().err


def test_gaff1_base_autodetect(tmp_path_factory, tmp_path):
    d = str(tmp_path_factory.mktemp("gaff1"))
    parametrize(ACETAMINOPHEN, d, atom_type="gaff")
    parm, rst = build(d, preamble="source leaprc.gaff")
    rc, report = run(parm, rst, tmp_path)
    assert rc == 0, report
    assert "base force field: gaff " in report


def test_ambiguous_autodetect(gaff2_box, tmp_path_factory, tmp_path, capsys):
    src, _, _ = gaff2_box
    d = str(tmp_path_factory.mktemp("two_ligs"))
    for f in ("LIG.mol2", "LIG.frcmod"):
        os.symlink(os.path.join(src, f), os.path.join(d, f))
    parm, rst = build(d, solute="L2 = copy LIG\ntranslate L2 {10 0 0}\nsys = combine {LIG L2}")
    rc, _ = run(parm, rst, tmp_path)
    assert rc == 2 and "ambiguous" in capsys.readouterr().err
    rc, report = run(parm, rst, tmp_path, "--mask", ":1", "--no-check")
    assert rc == 0, report


def test_multi_residue_mask_rejected(gaff2_box, tmp_path, capsys):
    _, parm, rst = gaff2_box
    rc, _ = run(parm, rst, tmp_path, "--mask", ":1-2")
    assert rc == 2 and "spans 2 residues" in capsys.readouterr().err


def test_nbfix_warns(gaff2_box, tmp_path_factory, tmp_path):
    from parmed.tools import changeLJPair
    _, parm, rst = gaff2_box
    p = pmd.load_file(parm, rst)
    changeLJPair(p, "@%oh", "@%OW", "3.5", "0.2").execute()
    nparm = str(tmp_path_factory.mktemp("nbfix") / "nbfix.parm7")
    p.save(nparm)
    rc, report = run(nparm, rst, tmp_path, "--no-check")
    assert rc == 0
    assert re.search(r"WARNING: .*off-diagonal", report)


def test_validation_catches_missing_terms(gaff2_box, tmp_path, monkeypatch):
    from unleap import frcmod
    real = frcmod.compare_to_base

    def drop_dihedrals(full, base):
        out = real(full, base)
        out["dihedral"] = []
        return out

    monkeypatch.setattr(frcmod, "compare_to_base", drop_dihedrals)
    _, parm, rst = gaff2_box
    rc, report = run(parm, rst, tmp_path)
    assert rc == 1
    assert "[FAIL] lib + diff frcmod" in report and "[PASS] lib + full frcmod" in report


def test_netcdf_trajectory_and_remd_rst7(gaff2_box, tmp_path_factory, tmp_path):
    import numpy as np
    from parmed.amber import NetCDFTraj
    _, parm, rst = gaff2_box
    d = tmp_path_factory.mktemp("coords")
    s = pmd.load_file(parm, rst)
    nc = NetCDFTraj.open_new(str(d / "traj.nc"), natom=len(s.atoms), box=True)
    for shift in (0.0, 0.5):
        nc.add_coordinates(np.array(s.coordinates) + shift)
        nc.add_box(s.box)
    nc.close()
    rc, report = run(parm, str(d / "traj.nc"), tmp_path / "nc", "--frame", "-1")
    assert rc == 0, report
    # pmemd/REMD restarts carry a temperature on the natom line, which ParmEd's sniffing rejects
    lines = open(rst).read().splitlines()
    lines[1] = lines[1].rstrip() + "  0.2980000E+03"
    with open(d / "remd.rst7", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    rc, report = run(parm, str(d / "remd.rst7"), tmp_path / "remd")
    assert rc == 0, report


def test_reversed_improper_order(gaff2_box):
    from unleap.frcmod import improper_atoms
    _, parm, rst = gaff2_box
    s = pmd.load_file(parm, rst)
    s.strip("!:LIG")
    for d in (d for d in s.dihedrals if d.improper):
        center = improper_atoms(d)[2]
        assert all(center in a.bond_partners for a in (d.atom1, d.atom2, d.atom3, d.atom4) if a is not center)
        d.atom1, d.atom2, d.atom3, d.atom4 = d.atom4, d.atom3, d.atom2, d.atom1   # tleap's reversed form
        assert improper_atoms(d)[2] is center


def test_leap_angle_correction_recovers_input(gaff2_box, tmp_path):
    d, parm, rst = gaff2_box
    rc, report = run(parm, rst, tmp_path)
    assert rc == 0, report
    # every equilibrium angle comes back with at most the 2-3 decimals GAFF2 uses
    for line in sections(tmp_path / "LIG.full.frcmod")["ANGL"]:
        theta = float(line.split()[-1])
        assert abs(theta - round(theta, 3)) < 1e-5, line


def test_sdf_is_rdkit_readable(charged_complex, tmp_path):
    from rdkit import Chem
    _, parm, rst = charged_complex
    rc, report = run(parm, rst, tmp_path, "--no-check")
    assert rc == 0, report
    mol = Chem.MolFromMolFile(str(tmp_path / "LIG.sdf"), removeHs=False)
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == Chem.CanonSmiles("[O-]C(=O)c1ccc(Cl)cc1")
