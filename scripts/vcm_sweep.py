"""Sweep unleap over ABFE/vcm_* jobs and compare against each job's inputs/abfe/<L>/binder_<L>.*

Usage: python scripts/vcm_sweep.py <ABFE dir containing vcm_*> <output dir>
"""
import glob
import json
import os
import subprocess
import sys

import parmed as pmd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from unleap import UnleapError                     # noqa: E402
from unleap.amber import run_tleap, sander_energy       # noqa: E402
from unleap.validate import compare_energies, compare_terms   # noqa: E402

ABFE, OUT = (os.path.abspath(os.path.expanduser(a)) for a in sys.argv[1:3])



def find_inputs(job, kind):
    """(parm7, coords) from the trial-1 production dir: its parm7 + the lambda=0 production rst7."""
    for d in sorted(glob.glob(f"{job}/results/abfe_run/abfe/attempt-*/*/t1_{kind}_abfeproduction"),
                    reverse=True):
        parms = glob.glob(f"{d}/{kind}_*.parm7")
        rst = sorted(r for r in glob.glob(f"{d}/{kind}_*_*.rst7"))   # excludes in_*/…ref_* restarts
        if parms and rst:
            return parms[0], rst[0]
    return None, None


def reference(job):
    mols = sorted(glob.glob(f"{job}/inputs/abfe/*/binder_*.mol2"))
    if not mols:
        mols = sorted(glob.glob(f"{job}/results/abfe_run/abfe/attempt-*/*/binder_*.mol2"))
    if not mols:
        return None
    base = mols[0][:-5]
    return {k: base + ext for k, ext in (("mol2", ".mol2"), ("frcmod", ".frcmod"), ("lib", ".lib"))
            if os.path.exists(base + ext)}


def lib_unit(path):
    lines = open(path).read().splitlines()
    return lines[lines.index("!!index array str") + 1].strip().strip('"')


def compare_ref(outdir, ref):
    """Rebuild the job's original ligand files and compare to what parm2lig extracted."""
    orig = glob.glob(f"{outdir}/*.orig.parm7")[0]
    orig_rst = orig[:-6] + ".rst7"
    work = os.path.join(outdir, "ref")
    os.makedirs(work, exist_ok=True)
    unit = f"m = loadmol2 {ref['mol2']}"
    if "lib" in ref:
        unit = f"loadoff {ref['lib']}\nm = copy {lib_unit(ref['lib'])}"
    frc = f"loadamberparams {ref['frcmod']}" if "frcmod" in ref else ""
    run_tleap(f"source leaprc.gaff2\n{frc}\n{unit}\nsaveamberparm m ref.parm7 ref.rst7", work, "ref")
    o = pmd.load_file(orig, orig_rst)
    r = pmd.load_file(os.path.join(work, "ref.parm7"))
    if [a.name for a in o.atoms] != [a.name for a in r.atoms]:
        return ["atom order/names differ from reference"]
    # masses: compare against the reference build (plain masses), i.e. checks the HMR undo too
    mine = pmd.load_file(glob.glob(f"{outdir}/check/rebuild1.parm7")[0])
    bad = compare_terms(o, r, [a.mass for a in mine.atoms])
    bad += compare_energies(sander_energy(orig, orig_rst, work, "o"),
                            sander_energy(os.path.join(work, "ref.parm7"), orig_rst, work, "r"))
    return bad


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for job in sorted(glob.glob(f"{ABFE}/vcm_*/")):
        job = job.rstrip("/")
        name = os.path.basename(job)
        ref = reference(job)
        for kind in ("complex", "binder"):
            parm, crd = find_inputs(job, kind)
            row = {"job": name, "kind": kind, "parm": parm and os.path.relpath(parm, job), "crd": crd and os.path.relpath(crd, job)}
            if parm is None:
                row["status"] = "no topology+coords"
                rows.append(row)
                continue
            outdir = os.path.join(OUT, f"{name}_{kind}")
            args = [sys.executable, "-m", "unleap.cli", "-p", parm, "-c", crd, "-o", outdir]
            if crd.endswith(".nc"):
                args += ["--frame", "-1"]
            res = subprocess.run(args, capture_output=True, text=True, cwd=REPO)
            report = open(f"{outdir}/report.txt").read() if os.path.exists(f"{outdir}/report.txt") else ""
            row["rc"] = res.returncode
            row["roundtrip"] = f"{report.count('[PASS]')}/3"
            row["hmr"] = "HMR detected" in report
            row["base"] = next((l.split()[3] for l in report.splitlines() if l.startswith("base force field")), None)
            row["warnings"] = [l for l in report.splitlines() if l.startswith("WARNING")]
            if res.returncode == 2:
                row["error"] = res.stderr.strip().splitlines()[-1]
            if res.returncode == 1:
                row["fail"] = [l.strip() for l in report.splitlines() if l.startswith("         ")][:6]
            if ref and res.returncode in (0, 1):
                try:
                    bad = compare_ref(outdir, ref)
                    row["vs_inputs"] = "identical" if not bad else bad[:6]
                except (UnleapError, Exception) as exc:
                    row["vs_inputs"] = f"ref build failed: {exc}"
            elif not ref:
                row["vs_inputs"] = "no reference inputs"
            rows.append(row)
            print(json.dumps(row), flush=True)
    json.dump(rows, open(os.path.join(OUT, "summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
