"""Recover a ligand's mol2 / frcmod / OFF lib from an Amber parm7 + rst7 pair."""

import logging

__version__ = "0.2.1"

__all__ = ["extract_ligand", "extract_ligand_from_structure", "LigandParams", "CheckResult",
           "BaseFF", "UnleapError", "__version__"]

logging.getLogger("unleap").addHandler(logging.NullHandler())


class UnleapError(RuntimeError):
    """Raised for user-facing failures (bad selection, unsupported topology, ...)."""


# the API pulls in parmed and rdkit; load it on first use so `unleap.__version__` stays cheap
_LAZY = {"extract_ligand": "api", "extract_ligand_from_structure": "api", "LigandParams": "api",
         "CheckResult": "validate", "BaseFF": "amber"}


def __getattr__(name):
    if name in _LAZY:
        import importlib
        return getattr(importlib.import_module(f".{_LAZY[name]}", __name__), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
