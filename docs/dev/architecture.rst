Architecture
============

Layout
------

.. code-block:: text

   unleap/
     __init__.py   version, UnleapError (user-facing failures -> exit code 2)
     cli.py        argparse + orchestration; writes report.txt
     amber.py      AMBERHOME helpers, base force-field presets, tleap / sander runners
     selection.py  --mask / --resname / auto-detect (standard residue name sets)
     extract.py    coordinate reading, ligand strip, sanity checks, PBC unwrap, HMR undo,
                   atom-type consistency
     frcmod.py     full parameter set, tleap-style base lookup, diff, base scoring, frcmod writer
     lib.py        OFF library via ParmEd ResidueTemplate
     mol2.py       Tripos mol2 + SDF writers, RDKit bond-order inference
     validate.py   tleap rebuilds, term-by-term comparison, sander energy comparison
   tests/
     conftest.py   builds fixtures on the fly (RDKit -> antechamber -c gas -> parmchk2 -> tleap)
     test_unleap.py
   scripts/
     vcm_sweep.py  real-data regression sweep over ABFE vcm_* jobs
   docs/           this documentation (Sphinx, Read the Docs theme)

Data flow
---------

:func:`unleap.cli.run` runs these steps in order:

1. :func:`unleap.extract.load` loads the topology with ParmEd and reads coordinates with
   :func:`~unleap.extract.read_coordinates` (see :doc:`amber_notes`, item 8).
2. :func:`unleap.selection.ligand_mask`, then :func:`unleap.extract.strip`. This checks the
   selection (exactly one whole residue, not covalent), warns about NBFIX/12-6-4/CMAP, then
   **reloads** the topology and ``strip("!(mask)")``\ s it, then unwraps the ligand across PBC.
3. :func:`unleap.frcmod.full_params` builds the parameter set with
   ``AmberParameterSet.from_structure``, then applies the masses, improper dedupe and order,
   and the leap angle correction.
4. :func:`unleap.cli.pick_base` scores each preset (or uses ``--base``) and returns the
   :class:`~unleap.amber.BaseFF` and its loaded ``.dat``.
5. :func:`unleap.extract.undo_hmr` (needs the base for reference H masses), then
   :func:`~unleap.extract.check_type_consistency`, then ``full_params`` again so the masses are final.
6. :func:`unleap.frcmod.diff_params` compares each term with the base using tleap's lookup rules
   (:doc:`diff_rules`).
7. The frcmods, lib, mol2/sdf and ``tleap.in`` are written.
8. :func:`unleap.validate.run` saves ``LIG.orig.parm7``, does the three tleap rebuilds, then runs
   :func:`~unleap.validate.compare_terms` and :func:`~unleap.validate.compare_energies`.

State between steps
-------------------

:class:`unleap.extract.Ligand` carries the state between steps. ``ligand.masses`` is the
per-atom mass list to write. It differs from ``lig.atoms[i].mass`` only when HMR was undone.

.. warning::

   Never mutate parameter types attached to ``ligand.lig``. :func:`unleap.validate.save_stripped`
   writes that structure as the reference everything is compared against. That's why
   :func:`unleap.frcmod.full_params` replaces angle and dihedral types with corrected *copies*.
