"""unleap: strip a ligand out of a parm7/rst7 pair and regenerate its mol2, frcmod and lib."""

import argparse
import os
import sys

from . import UnleapError, __version__, extract, frcmod, lib, mol2, selection, validate
from .amber import preset_bases, resolve_base


def build_parser():
    p = argparse.ArgumentParser(prog="unleap", description=__doc__)
    p.add_argument("-p", "--parm", required=True, help="input parm7/prmtop")
    p.add_argument("-c", "--crd", required=True,
                   help="coordinates: rst7/inpcrd, NetCDF restart, or NetCDF trajectory (see --frame)")
    p.add_argument("--frame", type=int, default=0,
                   help="frame to take from a NetCDF trajectory, 0-based; negative counts from the end "
                        "(default: %(default)s)")
    sel = p.add_mutually_exclusive_group()
    sel.add_argument("--mask", help="Amber mask selecting the ligand residue (e.g. ':LIG' or ':5')")
    sel.add_argument("--resname", help="ligand residue name")
    p.add_argument("-o", "--outdir", default="unleap_out", help="output directory (default: %(default)s)")
    p.add_argument("--name", help="unit/residue name for the outputs (default: the original residue name)")
    p.add_argument("--base", default="auto",
                   help="base force field for the diff frcmod: auto, gaff, gaff2, gaff211 or a .dat path "
                        "(default: %(default)s)")
    p.add_argument("--keep-masses", action="store_true", help="do not undo hydrogen mass repartitioning")
    p.add_argument("--no-check", action="store_true", help="skip the tleap/sander round-trip validation")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


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
        raise UnleapError("no gaff*.dat found in AMBERHOME; pass --base path/to/file.dat")
    return best[0], best[1], scores


def run(args):
    report = []
    log = lambda msg="": (print(msg), report.append(msg))
    os.makedirs(args.outdir, exist_ok=True)

    struct = extract.load(args.parm, args.crd, args.frame)
    mask = selection.ligand_mask(struct, args.mask, args.resname)
    ligand = extract.strip(struct, mask, args.parm)
    name = args.name or ligand.resname
    if args.name:
        ligand.lig.residues[0].name = name
    log(f"unleap {__version__}")
    log(f"input: {args.parm} / {args.crd}")
    log(f"ligand: mask {mask} -> residue {ligand.resname}, {len(ligand.lig.atoms)} atoms, "
        f"net charge {sum(a.charge for a in ligand.lig.atoms):+.4f}")

    full = frcmod.full_params(ligand)
    base, base_params, scores = pick_base(args.base, full)
    log("base force field: " + base.name + "  (match: "
        + ", ".join(f"{k} {v:.1%}" for k, v in scores.items()) + ")")

    if not args.keep_masses and extract.undo_hmr(ligand, base_params):
        h = [(a.mass, m) for a, m in zip(ligand.lig.atoms, ligand.masses) if a.atomic_number == 1]
        log(f"HMR detected: H masses {max(x for x, _ in h):.3f} -> {max(y for _, y in h):.3f}, "
            "excess returned to the bonded heavy atoms")
    extract.check_type_consistency(ligand)
    full = frcmod.full_params(ligand)
    diff, counts = frcmod.diff_params(full, base_params)

    files = {k: os.path.abspath(os.path.join(args.outdir, f"{name}{ext}")) for k, ext in
             (("mol2", ".mol2"), ("frcmod", ".frcmod"), ("full", ".full.frcmod"), ("lib", ".lib"))}
    frcmod.write(diff, files["frcmod"], f"unleap: {name} terms missing from/different in {base.name} "
                                        f"(load after: {base.leap_preamble})")
    frcmod.write(full, files["full"], f"unleap: all {name} parameters (self-contained)")
    lib.write(ligand, files["lib"], name)
    sdf = os.path.abspath(os.path.join(args.outdir, f"{name}.sdf"))
    if mol2.write(ligand, files["mol2"], name, sdf):
        log(f"wrote {sdf}")
    else:
        ligand.warnings.append("RDKit could not assign bond orders; mol2 bonds are all single")
    with open(os.path.join(args.outdir, "tleap.in"), "w") as fh:
        fh.write(f"{base.leap_preamble}\nloadamberparams {name}.frcmod\nloadoff {name}.lib\n")
    log("diff frcmod: " + ", ".join(f"{v} {k}" for k, v in counts.items()))
    for path in files.values():
        log(f"wrote {path}")

    ok = True
    if not args.no_check:
        results, _ = validate.run(ligand, args.outdir, name, base.leap_preamble, files)
        log()
        log("round-trip validation (vs stripped original):")
        for r in results:
            log(f"  [{'PASS' if r.ok else 'FAIL'}] {r.name}")
            for d in r.details[:20]:
                log(f"         {d}")
            if len(r.details) > 20:
                log(f"         ... {len(r.details) - 20} more")
            ok &= r.ok
    for w in ligand.warnings:
        log(f"WARNING: {w}")

    with open(os.path.join(args.outdir, "report.txt"), "w") as fh:
        fh.write("\n".join(report) + "\n")
    return 0 if ok else 1


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except UnleapError as exc:
        print(f"unleap: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
