"""Python API: extract a ligand's parameters from a built system and get them back as data."""

import logging
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from . import UnleapError, __version__, extract, frcmod, lib, mol2, selection, validate
from .amber import BaseFF, preset_bases, resolve_base

log = logging.getLogger("unleap")


@dataclass
class LigandParams:
    """Everything unleap recovers for one ligand. File bodies are strings; :meth:`write` puts
    them on disk."""

    name: str                   # unit/residue name used in the outputs
    resname: str                # residue name in the input topology
    mask: str
    n_atoms: int
    net_charge: float
    frcmod: str                 # diff against ``base``
    full_frcmod: str            # self-contained
    lib: str
    mol2: str
    sdf: object                 # str, or None when RDKit could not assign bond orders
    tleap_in: str
    base: BaseFF
    base_scores: dict
    hmr_undone: bool
    diff_counts: dict
    checks: list = field(default_factory=list)      # [CheckResult]; empty unless validated
    validated: bool = False
    warnings: list = field(default_factory=list)
    report: str = ""
    structure: object = field(default=None, repr=False)   # stripped ligand AmberParm, input masses

    @property
    def ok(self):
        """True when every validation check passed (vacuously True when not validated)."""
        return all(c.ok for c in self.checks)

    def files(self):
        """{file name: contents} for every text output."""
        out = {f"{self.name}.frcmod": self.frcmod, f"{self.name}.full.frcmod": self.full_frcmod,
               f"{self.name}.lib": self.lib, f"{self.name}.mol2": self.mol2, "tleap.in": self.tleap_in}
        if self.sdf is not None:
            out[f"{self.name}.sdf"] = self.sdf
        return out

    def write(self, outdir):
        """Write the outputs to ``outdir`` (created if needed); returns {kind: Path}.

        Kinds: frcmod, full_frcmod, lib, mol2, sdf (if any), tleap_in, and orig_parm7/orig_rst7
        (the stripped original ligand) when the result was validated.
        """
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        kinds = {f"{self.name}.frcmod": "frcmod", f"{self.name}.full.frcmod": "full_frcmod",
                 f"{self.name}.lib": "lib", f"{self.name}.mol2": "mol2", f"{self.name}.sdf": "sdf",
                 "tleap.in": "tleap_in"}
        paths = {}
        for fname, text in self.files().items():
            path = (outdir / fname).resolve()
            path.write_text(text)
            paths[kinds[fname]] = path
        if self.validated and self.structure is not None:
            parm, rst = validate.save_stripped(self.structure, str(outdir / f"{self.name}.orig"))
            paths["orig_parm7"], paths["orig_rst7"] = Path(parm).resolve(), Path(rst).resolve()
        return paths

    def to_dict(self):
        """JSON-able metadata (no file bodies)."""
        return {
            "unleap_version": __version__,
            "name": self.name, "resname": self.resname, "mask": self.mask,
            "n_atoms": self.n_atoms, "net_charge": self.net_charge,
            "base": {"name": self.base.name, "dat": self.base.dat, "leap_preamble": self.base.leap_preamble},
            "base_scores": dict(self.base_scores),
            "hmr_undone": self.hmr_undone,
            "diff_counts": dict(self.diff_counts),
            "has_sdf": self.sdf is not None,
            "validated": self.validated,
            "ok": self.ok,
            "checks": [c.to_dict() for c in self.checks],
            "warnings": list(self.warnings),
        }


def pick_base(spec, full):
    """Return (BaseFF, params, {name: score})."""
    base = resolve_base(spec)
    if base is not None:
        params = base.load()
        return base, params, {base.name: frcmod.score_base(full, params)}
    scores, best = {}, None
    for cand in preset_bases():
        params = cand.load()
        scores[cand.name] = frcmod.score_base(full, params)
        if best is None or scores[cand.name] > scores[best[0].name] + 1e-12:
            best = (cand, params)
    if best is None:
        raise UnleapError("no gaff*.dat found in AMBERHOME; pass base='path/to/file.dat'")
    return best[0], best[1], scores


def extract_ligand(parm7, crd, *, mask=None, resname=None, frame=0, name=None, base="auto",
                   undo_hmr=True, validate=True, workdir=None):
    """Recover the mol2 / frcmod / OFF lib of one ligand residue in a built Amber system.

    Parameters
    ----------
    parm7 : str or path
        Amber topology.
    crd : str or path
        Coordinates: ASCII rst7/inpcrd, NetCDF restart, or NetCDF trajectory (see ``frame``).
    mask, resname : str, optional
        Amber mask or residue name selecting the ligand (at most one). If both are None, the
        single non-standard residue is auto-detected.
    frame : int
        Trajectory frame, 0-based; negative counts from the end.
    name : str, optional
        Unit/residue name for the outputs (default: the original residue name).
    base : str
        Base force field for the diff frcmod: "auto", "gaff", "gaff2", "gaff211" or a .dat path.
    undo_hmr : bool
        Detect and reverse hydrogen mass repartitioning.
    validate : bool
        Rebuild the ligand with tleap and compare terms and sander energies to the original.
    workdir : str or path, optional
        Scratch directory for the validation (kept). A temporary one is used if None.

    Returns
    -------
    LigandParams

    Raises
    ------
    UnleapError
        For any user-facing failure (bad selection, unsupported topology, missing AmberTools).
    """
    struct = extract.load(os.fspath(parm7), os.fspath(crd), frame)
    return _extract(struct, f"{os.fspath(parm7)} / {os.fspath(crd)}", mask=mask, resname=resname,
                    name=name, base=base, undo_hmr=undo_hmr, validate=validate, workdir=workdir)


def extract_ligand_from_structure(struct, *, mask=None, resname=None, name=None, base="auto",
                                  undo_hmr=True, validate=True, workdir=None):
    """Like :func:`extract_ligand`, for an already-loaded ``parmed.amber.AmberParm`` with
    coordinates (and box, if periodic). ``struct`` is not modified, so it can be reused, e.g.
    once per resname in a dual-topology system. After editing it through ParmEd's object
    model, call ``struct.remake_parm()`` first."""
    extract.check_parm(struct)
    return _extract(struct, getattr(struct, "name", None) or "<structure>", mask=mask,
                    resname=resname, name=name, base=base, undo_hmr=undo_hmr, validate=validate,
                    workdir=workdir)


def _extract(struct, source, *, mask, resname, name, base, undo_hmr, validate, workdir):
    report = []
    say = lambda msg="": (log.info(msg), report.append(msg))

    mask = selection.ligand_mask(struct, mask, resname)
    ligand = extract.strip(struct, mask)
    orig_resname = ligand.resname
    name = name or orig_resname
    ligand.lig.residues[0].name = name
    net_charge = sum(a.charge for a in ligand.lig.atoms)
    say(f"unleap {__version__}")
    say(f"input: {source}")
    say(f"ligand: mask {mask} -> residue {orig_resname}, {len(ligand.lig.atoms)} atoms, "
        f"net charge {net_charge:+.4f}")

    full = frcmod.full_params(ligand)
    base_ff, base_params, scores = pick_base(base, full)
    say("base force field: " + base_ff.name + "  (match: "
        + ", ".join(f"{k} {v:.1%}" for k, v in scores.items()) + ")")

    if undo_hmr and extract.undo_hmr(ligand, base_params):
        h = [(a.mass, m) for a, m in zip(ligand.lig.atoms, ligand.masses) if a.atomic_number == 1]
        say(f"HMR detected: H masses {max(x for x, _ in h):.3f} -> {max(y for _, y in h):.3f}, "
            "excess returned to the bonded heavy atoms")
    extract.check_type_consistency(ligand)
    full = frcmod.full_params(ligand)
    diff, counts = frcmod.diff_params(full, base_params)

    texts = {
        "frcmod": frcmod.render(diff, f"unleap: {name} terms missing from/different in {base_ff.name} "
                                      f"(load after: {base_ff.leap_preamble})"),
        "full": frcmod.render(full, f"unleap: all {name} parameters (self-contained)"),
        "lib": lib.render(ligand, name),
    }
    texts["mol2"], sdf = mol2.render(ligand, name)
    if sdf is None:
        ligand.warnings.append("RDKit could not assign bond orders; mol2 bonds are all single")
    say("diff frcmod: " + ", ".join(f"{v} {k}" for k, v in counts.items()))

    checks = []
    if validate:
        if workdir is None:
            with tempfile.TemporaryDirectory(prefix="unleap_") as tmp:
                checks, _ = _validate(ligand, name, base_ff, texts, tmp)
        else:
            checks, _ = _validate(ligand, name, base_ff, texts, os.fspath(workdir))
        say()
        say("round-trip validation (vs stripped original):")
        for r in checks:
            say(f"  [{'PASS' if r.ok else 'FAIL'}] {r.name}")
            for d in r.details[:20]:
                say(f"         {d}")
            if len(r.details) > 20:
                say(f"         ... {len(r.details) - 20} more")
    for w in ligand.warnings:
        say(f"WARNING: {w}")

    return LigandParams(
        name=name, resname=orig_resname, mask=mask, n_atoms=len(ligand.lig.atoms), net_charge=net_charge,
        frcmod=texts["frcmod"], full_frcmod=texts["full"], lib=texts["lib"], mol2=texts["mol2"], sdf=sdf,
        tleap_in=f"{base_ff.leap_preamble}\nloadamberparams {name}.frcmod\nloadoff {name}.lib\n",
        base=base_ff, base_scores=scores, hmr_undone=ligand.hmr, diff_counts=counts,
        checks=checks, validated=bool(validate), warnings=list(ligand.warnings),
        report="\n".join(report) + "\n", structure=ligand.lig)


def _validate(ligand, name, base_ff, texts, workdir):
    return validate.run(ligand, name, base_ff.leap_preamble, texts, workdir)
