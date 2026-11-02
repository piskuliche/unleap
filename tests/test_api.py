import json
import os
import subprocess
import sys

import parmed as pmd
import pytest

from conftest import build
import unleap
from unleap.cli import main


def test_extract_ligand_matches_cli(gaff2_box, tmp_path, capsys):
    _, parm, rst = gaff2_box
    res = unleap.extract_ligand(parm, rst, workdir=tmp_path / "scratch")
    assert capsys.readouterr().out == ""        # the library never prints
    assert res.ok and res.validated and len(res.checks) == 3
    assert res.base.name == "gaff2" and res.resname == res.name == "LIG"
    assert res.sdf is not None and res.diff_counts["improper"] >= 0
    assert "[PASS] lib + full frcmod" in res.report

    paths = res.write(tmp_path / "api")
    assert main(["-p", parm, "-c", rst, "-o", str(tmp_path / "cli")]) == 0
    for kind, attr in (("frcmod", "frcmod"), ("full_frcmod", "full_frcmod"), ("lib", "lib"),
                       ("mol2", "mol2"), ("sdf", "sdf"), ("tleap_in", "tleap_in")):
        assert paths[kind].read_text() == getattr(res, attr)
        assert (tmp_path / "cli" / paths[kind].name).read_text() == getattr(res, attr), kind
    assert paths["orig_parm7"].exists() and paths["orig_rst7"].exists()


def test_to_dict_is_json_without_bodies(gaff2_box):
    _, parm, rst = gaff2_box
    res = unleap.extract_ligand(parm, rst)
    d = json.loads(json.dumps(res.to_dict()))
    assert d["unleap_version"] == unleap.__version__ and d["ok"] is True
    assert d["base"]["name"] == "gaff2" and [c["ok"] for c in d["checks"]] == [True] * 3
    assert "MASS" not in json.dumps(d)       # no frcmod text leaked in


def test_validate_false_needs_no_tleap(gaff2_box, monkeypatch):
    from unleap import amber, UnleapError

    def no_exe(name):
        raise UnleapError(f"{name} unavailable")

    monkeypatch.setattr(amber, "_exe", no_exe)
    _, parm, rst = gaff2_box
    res = unleap.extract_ligand(parm, rst, validate=False)
    assert res.checks == [] and not res.validated and res.ok
    assert "round-trip" not in res.report


def test_scratch_stays_in_workdir(gaff2_box, tmp_path, monkeypatch):
    _, parm, rst = gaff2_box
    cwd, work = tmp_path / "cwd", tmp_path / "work"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    unleap.extract_ligand(parm, rst, workdir=work)
    assert list(cwd.iterdir()) == []
    assert (work / "rebuild0.parm7").exists() and (work / "LIG.orig.parm7").exists()
    unleap.extract_ligand(parm, rst)                  # temporary workdir
    assert list(cwd.iterdir()) == []


def test_from_structure_leaves_input_untouched(charged_complex):
    _, parm, rst = charged_complex
    struct = pmd.load_file(parm, rst)
    natom, lj = len(struct.atoms), list(struct.parm_data["LENNARD_JONES_ACOEF"])
    res = unleap.extract_ligand_from_structure(struct, resname="LIG")
    assert len(struct.atoms) == natom and struct.parm_data["LENNARD_JONES_ACOEF"] == lj
    ref = unleap.extract_ligand(parm, rst)
    assert res.ok and (res.frcmod, res.lib, res.mol2) == (ref.frcmod, ref.lib, ref.mol2)
    again = unleap.extract_ligand_from_structure(struct, resname="LIG", validate=False)
    assert again.frcmod == res.frcmod


def test_from_structure_rejects_non_amber():
    with pytest.raises(unleap.UnleapError, match="AmberParm"):
        unleap.extract_ligand_from_structure(pmd.Structure())


def test_errors_raise_not_exit(gaff2_box, tmp_path_factory):
    src, _, _ = gaff2_box
    d = str(tmp_path_factory.mktemp("two_ligs_api"))
    for f in ("LIG.mol2", "LIG.frcmod"):
        os.symlink(os.path.join(src, f), os.path.join(d, f))
    parm, rst = build(d, solute="L2 = copy LIG\ntranslate L2 {10 0 0}\nsys = combine {LIG L2}")
    with pytest.raises(unleap.UnleapError, match="ambiguous"):
        unleap.extract_ligand(parm, rst)
    with pytest.raises(unleap.UnleapError, match="spans 2 residues"):
        unleap.extract_ligand(parm, rst, mask=":1-2")
    assert unleap.extract_ligand(parm, rst, mask=":2", name="L2X", validate=False).name == "L2X"


def test_version_import_is_light():
    code = "import sys, unleap; unleap.__version__; print('parmed' in sys.modules, 'rdkit' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.split() == ["False", "False"]
