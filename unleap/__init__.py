"""Recover a ligand's mol2 / frcmod / OFF lib from an Amber parm7 + rst7 pair."""

__version__ = "0.1.0"


class UnleapError(RuntimeError):
    """Raised for user-facing failures (bad selection, unsupported topology, ...)."""
