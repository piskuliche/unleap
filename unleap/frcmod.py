"""Build the full (self-contained) and diff-vs-base frcmod parameter sets."""

import math
from copy import copy
from itertools import permutations

import parmed as pmd

_TOL = dict(rel_tol=1e-5, abs_tol=1e-5)

# leap converts degrees with DEGTORAD = 0.0174533 (leap/basics.h), not pi/180, so every angle
# it writes is inflated by ~4.3e-7 when read back exactly; undo it so values round-trip.
LEAP_DEG_CORRECTION = (math.pi / 180.0) / 0.0174533


def improper_atoms(d):
    """Improper atoms ordered with the central atom third (frcmod convention).

    tleap may store an improper reversed in the prmtop (central atom second); a reversed
    torsion is the same angle, so reverse it back. The central atom is the one bonded to the
    other three.
    """
    atoms = (d.atom1, d.atom2, d.atom3, d.atom4)
    for n, center in enumerate(atoms):
        others = [a for a in atoms if a is not center]
        if all(center in a.bond_partners for a in others):
            if n == 1:
                return atoms[::-1]
            if n != 2:
                return tuple(others[:2]) + (center, others[2])
            return atoms
    return atoms


def _leap_inflated(values):
    """True if the angles look like leap's DEGTORAD output (rounder after correction)."""
    err = lambda c: sum(abs(v * c - round(v * c, 4)) for v in values)
    return bool(values) and err(LEAP_DEG_CORRECTION) < err(1.0)


def full_params(ligand):
    """Every parameter the ligand uses, with masses taken from ligand.masses (HMR undone)."""
    pset = pmd.amber.AmberParameterSet.from_structure(ligand.lig)
    for atom, mass in zip(ligand.lig.atoms, ligand.masses):
        pset.atom_types[atom.type].mass = mass
    # from_structure keys impropers under every outer-atom permutation; keep one, central third
    impropers = {}
    for d in ligand.lig.dihedrals:
        if d.improper:
            key = tuple(a.type for a in improper_atoms(d))
            impropers.setdefault(_improper_canon(key), (key, d.type))
    pset.improper_periodic_types = dict(impropers.values())

    corr = LEAP_DEG_CORRECTION if _leap_inflated([t.theteq for t in pset.angle_types.values()]) else 1.0
    ligand.angle_correction = corr
    if corr != 1.0:
        _rescale(pset, corr)
    return pset


def _rescale(pset, corr):
    """Replace angle/phase values with leap-corrected copies (never mutate the structure's types)."""
    new = {}
    for key, t in pset.angle_types.items():
        new.setdefault(id(t), pmd.AngleType(t.k, t.theteq * corr))
        pset.angle_types[key] = new[id(t)]
    fix = lambda t: pmd.DihedralType(t.phi_k, t.per, t.phase * corr, t.scee, t.scnb)
    for key, tl in pset.dihedral_types.items():
        if id(tl) not in new:
            lst = pmd.DihedralTypeList()
            for t in (tl if not isinstance(tl, pmd.DihedralType) else [tl]):
                lst.append(fix(t))
            new[id(tl)] = lst
        pset.dihedral_types[key] = new[id(tl)]
    for key, t in pset.improper_periodic_types.items():
        pset.improper_periodic_types[key] = fix(t)


def _improper_canon(key):
    return (key[2], tuple(sorted((key[0], key[1], key[3]))))


def _unique(table):
    """Iterate (key, value) once per term, skipping the reversed duplicate keys."""
    seen = set()
    for key, val in table.items():
        canon = min(key, key[::-1])
        if canon in seen:
            continue
        seen.add(canon)
        yield key, val


def _get(table, key):
    return table.get(key) or table.get(key[::-1])


def lookup_dihedral(base, key):
    """tleap order: explicit types (either direction), then X-b-c-X."""
    a, b, c, d = key
    return _get(base.dihedral_types, key) or _get(base.dihedral_types, ("X", b, c, "X"))


def lookup_improper(base, key):
    """Central atom is third; outer atoms may be permuted; most specific match wins."""
    a, b, c, d = key
    outer = (a, b, d)
    table = base.improper_periodic_types
    for nwild in range(3):
        for p in permutations(outer):
            cand = ("X",) * nwild + p[nwild:]
            if len(cand) == 3 and (cand[0], cand[1], c, cand[2]) in table:
                return table[(cand[0], cand[1], c, cand[2])]
    return None


def _close(x, y):
    return math.isclose(float(x), float(y), **_TOL)


def _dihedral_terms(dlist):
    if isinstance(dlist, pmd.DihedralType):
        dlist = [dlist]
    terms = []
    for t in dlist:
        if abs(t.phi_k) < 1e-8:
            continue
        terms.append((int(round(abs(t.per))), float(t.phase) % 360.0, float(t.phi_k),
                      float(t.scee), float(t.scnb)))
    return sorted(terms)


def same_dihedral(x, y, check_14=True):
    tx, ty = _dihedral_terms(x), _dihedral_terms(y)
    if len(tx) != len(ty):
        return False
    n = 5 if check_14 else 3
    for a, b in zip(tx, ty):
        if a[0] != b[0]:
            return False
        phase_diff = abs(a[1] - b[1]) % 360.0
        if min(phase_diff, 360.0 - phase_diff) > 1e-3:
            return False
        if not all(_close(u, v) for u, v in zip(a[2:n], b[2:n])):
            return False
    return True


def _same_type(lt, bt):
    return (bt is not None and _close(lt.mass, bt.mass) and _close(lt.rmin, bt.rmin)
            and _close(lt.epsilon, bt.epsilon))


def compare_to_base(full, base):
    """Return {section: [(key, value), ...]} of ligand terms that the base does not reproduce."""
    out = {"atom": [], "bond": [], "angle": [], "dihedral": [], "improper": []}
    for name, lt in full.atom_types.items():
        if not _same_type(lt, base.atom_types.get(name)):
            out["atom"].append((name, lt))
    for key, val in _unique(full.bond_types):
        bt = _get(base.bond_types, key)
        if bt is None or not (_close(val.k, bt.k) and _close(val.req, bt.req)):
            out["bond"].append((key, val))
    for key, val in _unique(full.angle_types):
        bt = _get(base.angle_types, key)
        if bt is None or not (_close(val.k, bt.k) and _close(val.theteq, bt.theteq)):
            out["angle"].append((key, val))
    for key, val in _unique(full.dihedral_types):
        bt = lookup_dihedral(base, key)
        if bt is None or not same_dihedral(val, bt):
            out["dihedral"].append((key, val))
    for key, val in full.improper_periodic_types.items():
        bt = lookup_improper(base, key)
        if bt is None or not same_dihedral(val, bt, check_14=False):
            out["improper"].append((key, val))
    return out


def count_terms(full):
    return {
        "atom": len(full.atom_types),
        "bond": sum(1 for _ in _unique(full.bond_types)),
        "angle": sum(1 for _ in _unique(full.angle_types)),
        "dihedral": sum(1 for _ in _unique(full.dihedral_types)),
        "improper": len(full.improper_periodic_types),
    }


def score_base(full, base):
    """Fraction of the ligand's bonded terms the base reproduces (masses ignored: HMR may alter them)."""
    diff = compare_to_base(full, base)
    total = count_terms(full)
    sections = ("bond", "angle", "dihedral", "improper")
    n = sum(total[s] for s in sections)
    bad = sum(len(diff[s]) for s in sections)
    lj_bad = sum(1 for name, lt in full.atom_types.items()
                 if base.atom_types.get(name) is None
                 or not (_close(lt.rmin, base.atom_types[name].rmin)
                         and _close(lt.epsilon, base.atom_types[name].epsilon)))
    return 1.0 - (bad + lj_bad) / max(1, n + total["atom"])


def diff_params(full, base):
    """Parameter set holding only the terms missing from, or different in, the base."""
    diff = compare_to_base(full, base)
    pset = pmd.amber.AmberParameterSet()
    for name, val in diff["atom"]:
        pset.atom_types[name] = copy(val)
    for section, table in (("bond", pset.bond_types), ("angle", pset.angle_types),
                           ("dihedral", pset.dihedral_types)):
        for key, val in diff[section]:
            table[key] = val
            table[key[::-1]] = val
    # All impropers are written, like parmchk2 does: when tleap matches a wildcard improper
    # instead of an explicit one it can store the four atoms in a different order, which changes
    # the improper angle (and energy) even though the parameters are identical.
    pset.improper_periodic_types = dict(full.improper_periodic_types)
    counts = {k: len(v) for k, v in diff.items()}
    counts["improper"] = len(full.improper_periodic_types)
    return pset, counts


def write(pset, path, title):
    """frcmod writer modeled on ParmEd's, but keeping prmtop precision (ParmEd rounds req/theteq
    to 3 decimals, which breaks GAFF2's 4-decimal bond lengths) and the prmtop improper order."""
    lines = [title.rstrip("\r\n"), "MASS"]
    for name, t in pset.atom_types.items():
        lines.append(f"{name:<6}{t.mass:8.4f}")
    lines += ["", "BOND"]
    for (a1, a2), t in _unique(pset.bond_types):
        lines.append(f"{a1:<2}-{a2:<2}   {t.k:12.6f}  {t.req:10.7f}")
    lines += ["", "ANGLE"]
    for (a1, a2, a3), t in _unique(pset.angle_types):
        lines.append(f"{a1:<2}-{a2:<2}-{a3:<2}   {t.k:12.6f}  {t.theteq:12.6f}")
    lines += ["", "DIHE"]
    for (a1, a2, a3, a4), tl in _unique(pset.dihedral_types):
        terms = [tl] if isinstance(tl, pmd.DihedralType) else list(tl)
        for n, t in enumerate(terms):
            per = abs(t.per) * (1 if n == len(terms) - 1 else -1)
            lines.append(f"{a1:<2}-{a2:<2}-{a3:<2}-{a4:<2} {1:4d} {t.phi_k:14.8f} {t.phase:8.3f} {per:5.1f}"
                         f"    SCEE={t.scee:g} SCNB={t.scnb:g}")
    lines += ["", "IMPROPER"]
    for (a1, a2, a3, a4), t in pset.improper_periodic_types.items():
        lines.append(f"{a1:<2}-{a2:<2}-{a3:<2}-{a4:<2} {t.phi_k:14.8f} {t.phase:8.3f} {t.per:5.1f}")
    lines += ["", "NONB"]
    for name, t in pset.atom_types.items():
        lines.append(f"{name:<2}  {t.rmin:12.8f} {t.epsilon:12.8f}")
    lines += ["", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
