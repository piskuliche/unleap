"""Sphinx configuration for unleap."""

import importlib
import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "unleap"
author = "Zeke Piskulich"
copyright = "2026, Zeke Piskulich"

from unleap import __version__ as release  # noqa: E402  (unleap/__init__.py has no heavy imports)

version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

def _missing(module):
    try:
        importlib.import_module(module)
        return False
    except ImportError:
        return True


# Read the Docs (and any machine without AmberTools) has no ParmEd/RDKit; mock them for autodoc.
autodoc_mock_imports = [m for m in ("parmed", "rdkit", "numpy", "scipy") if _missing(m)]
autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "undoc-members": False, "show-inheritance": True}
napoleon_numpy_docstring = True
napoleon_google_docstring = False

intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

templates_path = []
exclude_patterns = ["_build"]

html_theme = "sphinx_rtd_theme"
html_theme_options = {"navigation_depth": 3, "collapse_navigation": False}
html_title = f"unleap {release}"
