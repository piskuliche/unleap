"""Round-trip check: rebuild the ligand with tleap from the outputs and compare to the original."""

import copy
import logging
import math
import os
from dataclasses import dataclass, field

import parmed as pmd

from . import UnleapError
from .amber import run_tleap, sander_energy
from .frcmod import improper_atoms, same_dihedral

log = logging.getLogger("unleap")

ENERGY_TERMS = ("BOND", "ANGLE", "DIHED", "VDWAALS", "EEL", "1-4 VDW", "1-4 EEL")


@dataclass
class CheckResult:
    """One rebuild-and-compare check: ``details`` lists each mismatch (empty when ``ok``)."""

    name: str
    ok: bool
    details: list = field(default_factory=list)

    def to_dict(self):
        return {"name": self.name, "ok": self.ok, "details": list(self.details)}


def save_stripped(struct, prefix):
    """Write the stripped original ligand (non-periodic) for comparison."""
    lig = copy.copy(struct)
    lig.box = None
    lig.save(prefix + ".parm7", overwrite=True)
    lig.save(prefix + ".rst7", format="rst7", overwrite=True)
    return prefix + ".parm7", prefix + ".rst7"


def _terms(struct):
    t = {"bond": {}, "angle": {}, "dihedral": {}, "improper": {}}
    for b in struct.bonds:
        if b.type is not None:
            t["bond"][tuple(sorted((b.atom1.idx, b.atom2.idx)))] = (b.type.k, b.type.req)
    for a in struct.angles:
        i, j, k = a.atom1.idx, a.atom2.idx, a.atom3.idx
        t["angle"][(min(i, k), j, max(i, k))] = (a.type.k, a.type.theteq)
    for d in struct.dihedrals:
        q = (d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)
        if d.improper:
            c = [a.idx for a in improper_atoms(d)]
            key = (c[2], tuple(sorted((c[0], c[1], c[3]))))
            t["improper"].setdefault(key, []).append(d.type)
        else:
            t["dihedral"].setdefault(min(q, q[::-1]), []).append(d.type)
    return t


def _close(x, y, tol=1e-5):
    return math.isclose(x, y, rel_tol=tol, abs_tol=tol)


def compare_terms(orig, rebuilt, masses):
    """Term-by-term comparison keyed by atom indices. Returns a list of mismatch strings."""
    bad = []
    if len(orig.atoms) != len(rebuilt.atoms):
        return [f"atom count {len(orig.atoms)} != {len(rebuilt.atoms)}"]
    for a, b, m in zip(orig.atoms, rebuilt.atoms, masses):
        for what, x, y in (("charge", a.charge, b.charge), ("rmin", a.rmin, b.rmin),
                           ("epsilon", a.epsilon, b.epsilon), ("mass", m, b.mass)):
            if not _close(x, y, 1e-4):
                bad.append(f"atom {a.name}: {what} {x:.6f} != {y:.6f}")
        if a.name != b.name or a.type != b.type:
            bad.append(f"atom {a.idx + 1}: {a.name}/{a.type} != {b.name}/{b.type}")
    to, tr = _terms(orig), _terms(rebuilt)
    names = [a.name for a in orig.atoms]
    label = lambda idx: "-".join(names[i] for i in idx)
    for sec in ("bond", "angle"):
        for key in to[sec].keys() | tr[sec].keys():
            x, y = to[sec].get(key), tr[sec].get(key)
            if x is None or y is None:
                bad.append(f"{sec} {label(key)}: {'missing' if y is None else 'extra'} in rebuild")
            elif not all(_close(u, v) for u, v in zip(x, y)):
                bad.append(f"{sec} {label(key)}: {x} != {y}")
    for sec in ("dihedral", "improper"):
        for key in to[sec].keys() | tr[sec].keys():
            x, y = to[sec].get(key), tr[sec].get(key)
            lab = label(key) if sec == "dihedral" else label((key[1][0], key[1][1], key[0], key[1][2]))
            if x is None or y is None:
                bad.append(f"{sec} {lab}: {'missing' if y is None else 'extra'} in rebuild")
            elif not same_dihedral(x, y, check_14=(sec == "dihedral")):
                bad.append(f"{sec} {lab}: terms differ")
    return bad


def compare_energies(e_orig, e_new):
    bad = []
    for term in ENERGY_TERMS:
        x, y = e_orig.get(term), e_new.get(term)
        if x is None or y is None:
            continue
        if abs(x - y) > 2e-3 + 1e-6 * abs(x):
            bad.append(f"{term}: {x:.4f} != {y:.4f}")
    return bad


def run(ligand, name, preamble, texts, workdir):
    """Rebuild three ways in ``workdir`` and compare each to the stripped original.

    ``texts`` maps "frcmod", "full", "lib" and "mol2" to the rendered file contents.
    Returns ([CheckResult], original energies).
    """
    os.makedirs(workdir, exist_ok=True)
    files = {}
    for key, ext in (("frcmod", ".frcmod"), ("full", ".full.frcmod"), ("lib", ".lib"), ("mol2", ".mol2")):
        files[key] = os.path.abspath(os.path.join(workdir, name + ext))
        with open(files[key], "w") as fh:
            fh.write(texts[key])
    orig_parm, orig_rst = save_stripped(ligand.lig, os.path.join(workdir, f"{name}.orig"))
    orig = pmd.load_file(orig_parm, orig_rst)
    e_orig = sander_energy(orig_parm, orig_rst, workdir, "orig")

    variants = {
        "lib + diff frcmod": f"{preamble}\nloadamberparams {files['frcmod']}\nloadoff {files['lib']}\nm = copy {name}",
        "lib + full frcmod": f"loadamberparams {files['full']}\nloadoff {files['lib']}\nm = copy {name}",
        "mol2 + diff frcmod": f"{preamble}\nloadamberparams {files['frcmod']}\nm = loadmol2 {files['mol2']}",
    }
    results = []
    for n, (label, body) in enumerate(variants.items()):
        log.info("rebuilding: %s", label)
        tag = f"rebuild{n}"
        parm, rst = os.path.join(workdir, f"{tag}.parm7"), os.path.join(workdir, f"{tag}.rst7")
        try:
            run_tleap(f"{body}\nsaveamberparm m {parm} {rst}", workdir, tag)
            rebuilt = pmd.load_file(parm, rst)
            details = compare_terms(orig, rebuilt, ligand.masses)
            # score the rebuilt topology on the original coordinates: the lib/mol2 store positions
            # at ~1e-4 A, which alone shifts stiff angle energies by ~1e-3 kcal/mol
            details += compare_energies(e_orig, sander_energy(parm, orig_rst, workdir, tag))
        except UnleapError as exc:
            details = [str(exc)]
        results.append(CheckResult(label, not details, details))
    return results, e_orig
