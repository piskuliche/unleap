"""AmberTools environment helpers: base force fields, tleap and sander runners."""

import os
import re
import shutil
import subprocess
from dataclasses import dataclass

from . import UnleapError


@dataclass
class BaseFF:
    """A base force field the diff frcmod is written against."""

    name: str
    dat: str            # path to the .dat file
    leap_preamble: str  # tleap lines that load it

    def load(self):
        import parmed as pmd
        return pmd.amber.AmberParameterSet(self.dat)


def amberhome():
    home = os.environ.get("AMBERHOME")
    if not home or not os.path.isdir(home):
        raise UnleapError("AMBERHOME is not set; source AmberTools first (e.g. `source $AMBERHOME/amber.sh`)")
    return home


def parm_dir():
    return os.path.join(amberhome(), "dat", "leap", "parm")


# name -> (dat file, tleap preamble). gaff211 has no leaprc of its own.
_PRESETS = {
    "gaff": ("gaff.dat", "source leaprc.gaff"),
    "gaff2": ("gaff2.dat", "source leaprc.gaff2"),
    "gaff211": ("gaff211.dat", "loadamberparams gaff211.dat"),
}


def preset_bases():
    """All preset base force fields available in this AMBERHOME."""
    out = []
    for name, (dat, pre) in _PRESETS.items():
        path = os.path.join(parm_dir(), dat)
        if os.path.isfile(path):
            out.append(BaseFF(name, path, pre))
    return out


def resolve_base(spec):
    """Resolve a --base value (preset name or .dat path); 'auto' returns None."""
    if spec == "auto":
        return None
    if spec in _PRESETS:
        for b in preset_bases():
            if b.name == spec:
                return b
        raise UnleapError(f"{_PRESETS[spec][0]} not found in {parm_dir()}")
    if os.path.isfile(spec):
        path = os.path.abspath(spec)
        return BaseFF(os.path.splitext(os.path.basename(path))[0], path, f"loadamberparams {path}")
    raise UnleapError(f"unknown base force field {spec!r} (use auto, gaff, gaff2, gaff211 or a .dat path)")


def _exe(name):
    path = shutil.which(name) or os.path.join(amberhome(), "bin", name)
    if not os.path.isfile(path):
        raise UnleapError(f"{name} not found; source AmberTools first (e.g. `source $AMBERHOME/amber.sh`)")
    return path


def run_tleap(script, workdir, tag="tleap"):
    """Run a tleap script in workdir; return the log text. Raises on leap errors."""
    inp = os.path.join(workdir, f"{tag}.in")
    with open(inp, "w") as fh:
        fh.write(script.rstrip() + "\nquit\n")
    res = subprocess.run([_exe("tleap"), "-f", f"{tag}.in"], cwd=workdir,
                         capture_output=True, text=True)
    log = res.stdout + res.stderr
    with open(os.path.join(workdir, f"{tag}.log"), "w") as fh:
        fh.write(log)
    m = re.search(r"Errors = (\d+)", log)
    if res.returncode != 0 or m is None or int(m.group(1)) > 0 or "FATAL" in log:
        raise UnleapError(f"tleap failed; see {os.path.join(workdir, tag + '.log')}")
    return log


_MDIN = """single point energy
&cntrl
  imin=1, maxcyc=1, ncyc=1, ntb=0, igb=0, cut=999.0, ntpr=1, ntxo=1,
/
"""

_TERM_RE = re.compile(r"(1-4 VDW|1-4 EEL|1-4 NB|[A-Z][A-Z0-9]*)\s*=\s*(-?\d+\.\d+)")


def sander_energy(parm7, rst7, workdir, tag):
    """Single-point vacuum energy decomposition {term: kcal/mol} via sander."""
    mdin = os.path.join(workdir, f"{tag}.mdin")
    mdout = os.path.join(workdir, f"{tag}.mdout")
    with open(mdin, "w") as fh:
        fh.write(_MDIN)
    res = subprocess.run([_exe("sander"), "-O", "-i", mdin, "-o", mdout, "-p", parm7, "-c", rst7,
                          "-r", os.path.join(workdir, f"{tag}.ncrst"),
                          "-inf", os.path.join(workdir, f"{tag}.mdinfo")],
                         cwd=workdir, capture_output=True, text=True)
    if res.returncode != 0 or not os.path.isfile(mdout):
        raise UnleapError(f"sander failed for {parm7}: {res.stderr.strip()}")
    with open(mdout) as fh:
        text = fh.read()
    # first energy block after the NSTEP header is the energy at the input coordinates
    start = text.find("NSTEP       ENERGY")
    if start < 0:
        raise UnleapError(f"could not parse energies from {mdout}")
    block = text[start:].split("\n\n", 2)
    energies = {}
    for line in "\n".join(block[:2]).splitlines():
        for key, val in _TERM_RE.findall(line):
            energies.setdefault(key, float(val))
    energies.pop("NSTEP", None)
    return energies
