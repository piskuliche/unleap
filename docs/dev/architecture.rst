Architecture
============

Layout
------

.. code-block:: text

   unleap/
     __init__.py   version, UnleapError, lazy re-exports of the public API
     api.py        extract_ligand / extract_ligand_from_structure, LigandParams, base picking
     cli.py        argparse wrapper around api.extract_ligand; writes the files and report.txt
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
     test_unleap.py  CLI behaviour
     test_api.py     Python API: strings, to_dict, workdir, no printing, structure input
   scripts/
     vcm_sweep.py  real-data regression sweep over ABFE vcm_* jobs
   docs/           this documentation (Sphinx, Read the Docs theme)

Data flow
---------

:func:`unleap.extract_ligand` loads the files, then
:func:`unleap.extract_ligand_from_structure` runs the steps below. Both return a
:class:`~unleap.LigandParams`. The CLI (:func:`unleap.cli.run`) only calls ``extract_ligand``,
writes the result with :meth:`~unleap.LigandParams.write` and prints its report.

1. :func:`unleap.extract.load` loads the topology with ParmEd and reads coordinates with
   :func:`~unleap.extract.read_coordinates` (see :doc:`amber_notes`, item 8).
2. :func:`unleap.selection.ligand_mask`, then :func:`unleap.extract.strip`. This checks the
   selection (exactly one whole residue, not covalent), warns about NBFIX/12-6-4/CMAP, then
   ``strip("!(mask)")``\ s a ``copy.copy`` of the topology (AmberParm rebuilds the copy from
   ``parm_data``, so ``strip`` can rebuild the LJ tables), then unwraps the ligand across PBC.
3. :func:`unleap.frcmod.full_params` builds the parameter set with
   ``AmberParameterSet.from_structure``, then applies the masses, improper dedupe and order,
   and the leap angle correction.
4. :func:`unleap.api.pick_base` scores each preset (or uses ``--base``) and returns the
   :class:`~unleap.amber.BaseFF` and its loaded ``.dat``.
5. :func:`unleap.extract.undo_hmr` (needs the base for reference H masses), then
   :func:`~unleap.extract.check_type_consistency`, then ``full_params`` again so the masses are final.
6. :func:`unleap.frcmod.diff_params` compares each term with the base using tleap's lookup rules
   (:doc:`diff_rules`).
7. The frcmods, lib, mol2/sdf are rendered to strings (``frcmod.render``, ``lib.render``,
   ``mol2.render``); nothing touches the output directory.
8. :func:`unleap.validate.run` writes those strings and ``LIG.orig.parm7`` into the scratch
   ``workdir`` (a temporary directory unless one is given), does the three tleap rebuilds, then
   runs :func:`~unleap.validate.compare_terms` and :func:`~unleap.validate.compare_energies`.

The library reports progress through ``logging.getLogger("unleap")`` and never prints or exits;
the report text is returned as ``LigandParams.report``.

State between steps
-------------------

:class:`unleap.extract.Ligand` carries the state between steps. ``ligand.masses`` is the
per-atom mass list to write. It differs from ``lig.atoms[i].mass`` only when HMR was undone.

.. warning::

   Never mutate parameter types attached to ``ligand.lig``. :func:`unleap.validate.save_stripped`
   writes that structure as the reference everything is compared against. That's why
   :func:`unleap.frcmod.full_params` replaces angle and dihedral types with corrected *copies*.
