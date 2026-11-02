"""Load the system, strip out the ligand, and check it can be expressed as mol2/frcmod/lib."""

import math
from dataclasses import dataclass, field

import numpy as np
import parmed as pmd

from . import UnleapError

H_STANDARD_MASS = 1.008


@dataclass
class Ligand:
    struct: object              # full system (AmberParm)
    lig: object                 # stripped ligand (AmberParm), coordinates unwrapped
    mask: str
    hmr: bool = False
    masses: list = field(default_factory=list)   # per-atom masses to write (HMR undone if requested)
    warnings: list = field(default_factory=list)

    @property
    def resname(self):
        return self.lig.residues[0].name


def read_coordinates(path, frame=0):
    """(xyz, box) from an ASCII rst7/inpcrd or one frame of a NetCDF restart/trajectory.

    Reads the file directly: ParmEd's format sniffing rejects rst7 headers that carry a
    temperature (pmemd/REMD restarts), and loading a whole trajectory just for one frame is wasteful.
    """
    with open(path, "rb") as fh:
        is_netcdf = fh.read(3) == b"CDF"
    try:
        if not is_netcdf:
            rst = pmd.amber.Rst7.open(path)
            return np.asarray(rst.coordinates).reshape(-1, 3), (rst.box if rst.hasbox else None)
        nc = pmd.amber.NetCDFTraj.open_old(path) if _is_traj(path) else None
        if nc is None:
            rst = pmd.amber.NetCDFRestart.open_old(path)
            return np.asarray(rst.coordinates).reshape(-1, 3), (rst.box if rst.hasbox else None)
        var = nc._ncfile.variables
        nframes = var["coordinates"].shape[0]
        if not -nframes <= frame < nframes:
            raise UnleapError(f"--frame {frame} out of range; {path} has {nframes} frames")
        box = None
        if "cell_lengths" in var:
            box = list(np.asarray(var["cell_lengths"][frame])) + list(np.asarray(var["cell_angles"][frame]))
        return np.array(var["coordinates"][frame], dtype=float), box
    except UnleapError:
        raise
    except Exception as exc:
        raise UnleapError(f"could not read coordinates from {path}: {exc}") from exc


def _is_traj(path):
    """NetCDF trajectories have a 'frame' dimension; NetCDF restarts do not."""
    from scipy.io import netcdf_file
    with netcdf_file(path, "r", mmap=False) as nc:
        return "frame" in nc.dimensions


def load_parm(parm7):
    struct = pmd.load_file(parm7)
    if not isinstance(struct, pmd.amber.AmberParm) or isinstance(struct, pmd.amber.ChamberParm):
        raise UnleapError(f"{parm7} is not a plain Amber topology (Chamber/AMOEBA are unsupported)")
    return struct


def load(parm7, crd, frame=0):
    struct = load_parm(parm7)
    xyz, box = read_coordinates(crd, frame)
    if len(xyz) != len(struct.atoms):
        raise UnleapError(f"{crd} has {len(xyz)} atoms but {parm7} has {len(struct.atoms)}")
    struct.coordinates = xyz
    struct.box = box
    return struct


def strip(struct, mask, parm7):
    """Extract the ligand and run the topology sanity checks."""
    sel = pmd.amber.AmberMask(struct, mask).Selection()
    idx = {i for i, s in enumerate(sel) if s}
    if not idx:
        raise UnleapError(f"mask {mask!r} selects no atoms")
    residues = {struct.atoms[i].residue.idx for i in idx}
    if len(residues) != 1:
        raise UnleapError(f"mask {mask!r} spans {len(residues)} residues; only single-residue ligands are supported")
    res = struct.residues[residues.pop()]
    if len(res.atoms) != len(idx):
        raise UnleapError(f"mask {mask!r} selects part of residue {res.name}{res.number}; select the whole residue")
    for b in struct.bonds:
        if (b.atom1.idx in idx) != (b.atom2.idx in idx):
            raise UnleapError(f"ligand atom {b.atom1.name if b.atom1.idx in idx else b.atom2.name} is bonded "
                                   "outside the selection (covalent ligand); this is unsupported")

    # strip a fresh copy: AmberParm.strip rebuilds the LJ tables, Structure slicing does not
    lig = load_parm(parm7)
    lig.coordinates = struct.coordinates
    lig.box = struct.box
    lig.strip(f"!({mask})")
    out = Ligand(struct, lig, mask)
    out.warnings += _nonbonded_checks(struct, idx)
    if struct.cmaps and any(a.idx in idx for c in struct.cmaps for a in (c.atom1, c.atom3, c.atom5)):
        out.warnings.append("ligand has CMAP terms; these cannot be written to a frcmod and are dropped")
    if struct.box is not None:
        _unwrap(lig, struct.box)
    out.masses = [a.mass for a in lig.atoms]
    return out


def _nonbonded_checks(struct, idx):
    """Warn about LJ features a frcmod cannot express (NBFIX pairs, 12-6-4 C coefficients)."""
    warnings = []
    data = struct.parm_data
    ntypes = struct.pointers["NTYPES"]
    lig_types = sorted({data["ATOM_TYPE_INDEX"][i] for i in idx})
    radius, depth = struct.LJ_radius, struct.LJ_depth
    nbfix = set()
    for ti in lig_types:
        for tj in range(1, ntypes + 1):
            k = data["NONBONDED_PARM_INDEX"][ntypes * (ti - 1) + tj - 1]
            if k <= 0:
                continue
            eps = math.sqrt(depth[ti - 1] * depth[tj - 1])
            rmin = radius[ti - 1] + radius[tj - 1]
            acoef, bcoef = data["LENNARD_JONES_ACOEF"][k - 1], data["LENNARD_JONES_BCOEF"][k - 1]
            if not (math.isclose(acoef, eps * rmin ** 12, rel_tol=1e-4, abs_tol=1e-6)
                    and math.isclose(bcoef, 2 * eps * rmin ** 6, rel_tol=1e-4, abs_tol=1e-6)):
                nbfix.add(tj)
    if nbfix:
        warnings.append(f"ligand LJ types have {len(nbfix)} off-diagonal (NBFIX/changeLJPair) pair(s); "
                        "these are not written to the frcmod")
    ccoef = data.get("LENNARD_JONES_CCOEF")
    if ccoef:
        for ti in lig_types:
            for tj in range(1, ntypes + 1):
                k = data["NONBONDED_PARM_INDEX"][ntypes * (ti - 1) + tj - 1]
                if k > 0 and abs(ccoef[k - 1]) > 1e-12:
                    warnings.append("ligand has 12-6-4 C coefficients; these are not written to the frcmod")
                    return warnings
    return warnings


def _unwrap(lig, box):
    """Make bonded atoms contiguous under periodic boundaries (bond walk + minimum image)."""
    from parmed.geometry import box_lengths_and_angles_to_vectors
    vecs = np.array([[v.value_in_unit(pmd.unit.angstrom) if hasattr(v, "value_in_unit") else v
                      for v in vec] for vec in box_lengths_and_angles_to_vectors(*box)])
    inv = np.linalg.inv(vecs)
    xyz = np.array(lig.coordinates)
    done = [False] * len(lig.atoms)
    for root in range(len(lig.atoms)):
        if done[root]:
            continue
        done[root] = True
        stack = [root]
        while stack:
            i = stack.pop()
            for nb in lig.atoms[i].bond_partners:
                j = nb.idx
                if done[j]:
                    continue
                d = xyz[j] - xyz[i]
                frac = d @ inv
                xyz[j] = xyz[i] + (frac - np.round(frac)) @ vecs
                done[j] = True
                stack.append(j)
    lig.coordinates = xyz


def undo_hmr(ligand, base_params=None):
    """Detect hydrogen mass repartitioning and restore the original masses in ligand.masses.

    Each H gets its type's mass from the base force field (or 1.008); the mass it carried
    beyond that is returned to the heavy atom it is bonded to, which reverses HMR exactly.
    """
    atoms = ligand.lig.atoms
    if not any(a.atomic_number == 1 and a.mass > 1.5 for a in atoms):
        return False
    masses = [a.mass for a in atoms]
    for a in atoms:
        if a.atomic_number != 1:
            continue
        ref = H_STANDARD_MASS
        if base_params is not None and a.type in base_params.atom_types:
            ref = base_params.atom_types[a.type].mass or H_STANDARD_MASS
        excess = a.mass - ref
        masses[a.idx] = ref
        heavy = [p for p in a.bond_partners if p.atomic_number != 1]
        if len(heavy) == 1:
            masses[heavy[0].idx] += excess
        elif abs(excess) > 1e-6:
            ligand.warnings.append(f"cannot return HMR mass of {a.name}: not bonded to exactly one heavy atom")
    ligand.masses = [round(m, 4) for m in masses]
    ligand.hmr = True
    return True


def check_type_consistency(ligand):
    """Each atom type must map to one mass/LJ pair for a frcmod to describe it."""
    seen = {}
    for a, m in zip(ligand.lig.atoms, ligand.masses):
        key = (round(m, 4), round(a.rmin, 6), round(a.epsilon, 6))
        prev = seen.setdefault(a.type, key)
        if prev != key:
            hint = ""
            if prev[1:] == key[1:] and any(x.atomic_number == 1 and m > 1.5
                                           for x, m in zip(ligand.lig.atoms, ligand.masses)):
                hint = " (hydrogen mass repartitioning? drop --keep-masses to undo it)"
            raise UnleapError(f"atom type {a.type!r} has inconsistent mass/LJ within the ligand "
                                   f"({prev} vs {key}); it cannot be written as one frcmod type{hint}")
    return {t: k for t, k in seen.items()}
