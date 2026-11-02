Python API
==========

The command line is a thin wrapper around one function, :func:`unleap.extract_ligand`. Call it
from Python to get the parameters back as data (strings and plain values) instead of files.

.. code-block:: python

   import unleap

   res = unleap.extract_ligand("complex.parm7", "complex.rst7")   # ligand auto-detected
   res.ok              # True when all three round-trip checks passed
   res.frcmod          # diff frcmod text (against res.base)
   res.full_frcmod     # self-contained frcmod text
   res.lib, res.mol2   # OFF library and mol2 text
   res.sdf             # SDF text, or None when RDKit could not assign bond orders
   res.tleap_in        # tleap lines that load the outputs
   res.report          # the same text the CLI prints

   res.write("lig_out")   # what the CLI does: {kind: Path} of the files written

Arguments
---------

``extract_ligand(parm7, crd, *, mask=None, resname=None, frame=0, name=None, base="auto",
undo_hmr=True, validate=True, workdir=None)``

These match the command line options (:doc:`options`): ``undo_hmr=False`` is ``--keep-masses``
and ``validate=False`` is ``--no-check``. ``crd`` can be an rst7/inpcrd, a NetCDF restart or a
NetCDF trajectory (``frame`` picks the frame).

Results
-------

:class:`~unleap.LigandParams` also carries:

* ``name``, ``resname``, ``mask``, ``n_atoms``, ``net_charge``
* ``base`` (a :class:`~unleap.BaseFF`: ``name``, ``dat``, ``leap_preamble``) and ``base_scores``,
  the fraction of ligand terms each candidate base matched
* ``hmr_undone`` and ``diff_counts`` (number of terms per frcmod section)
* ``checks``: a list of :class:`~unleap.CheckResult` (``name``, ``ok``, ``details``)
* ``warnings``: things the frcmod cannot express (NBFIX, CMAP, ...), as strings
* ``structure``: the stripped ligand as a ParmEd ``AmberParm``

``res.to_dict()`` gives JSON-able metadata, including ``unleap_version``, without the file bodies.
That makes it a good provenance record to store next to the parameters.

Validation and scratch files
----------------------------

With ``validate=True`` the ligand is rebuilt with tleap and scored with sander
(:doc:`validation`). The scratch files go to ``workdir``, or to a temporary directory that is
deleted afterwards. Nothing is written to the current directory. ``validate=False`` skips the
rebuild, so tleap and sander are not needed. The base force field is still read from
``$AMBERHOME/dat/leap/parm``, unless ``base`` is a ``.dat`` path.

A failed check does not raise: inspect ``res.ok`` and ``res.checks``. Problems that stop the
extraction (bad selection, unsupported topology, missing AmberTools) raise
:class:`~unleap.UnleapError`.

Already-loaded structures
-------------------------

:func:`unleap.extract_ligand_from_structure` takes a ParmEd ``AmberParm`` that has coordinates.
The structure is not modified, so you can load a large system once and extract several ligands
from it. For example, both end states of a dual-topology RBFE system:

.. code-block:: python

   import parmed as pmd
   import unleap

   struct = pmd.load_file("unisc.parm7", "0.00000000_ti.rst7")
   results = {rn: unleap.extract_ligand_from_structure(struct, resname=rn) for rn in ("L00", "L02")}

If you edited the structure through ParmEd's object model, call ``struct.remake_parm()`` first:
unleap works from a copy that is rebuilt from the raw topology data.

Logging
-------

unleap never prints. Progress messages go to the ``unleap`` logger:

.. code-block:: python

   import logging
   logging.basicConfig(level=logging.INFO)
