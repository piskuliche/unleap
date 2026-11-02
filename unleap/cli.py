"""unleap: strip a ligand out of a parm7/rst7 pair and regenerate its mol2, frcmod and lib."""

import argparse
import os
import sys

from . import UnleapError, __version__
from .api import extract_ligand


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


def run(args):
    res = extract_ligand(args.parm, args.crd, mask=args.mask, resname=args.resname, frame=args.frame,
                         name=args.name, base=args.base, undo_hmr=not args.keep_masses,
                         validate=not args.no_check, workdir=os.path.join(args.outdir, "check"))
    paths = res.write(args.outdir)
    wrote = "".join(f"wrote {paths[k]}\n" for k in ("sdf", "mol2", "frcmod", "full_frcmod", "lib")
                    if k in paths)
    report = res.report + wrote
    print(report, end="")
    with open(os.path.join(args.outdir, "report.txt"), "w") as fh:
        fh.write(report)
    return 0 if res.ok else 1


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except UnleapError as exc:
        print(f"unleap: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
