Extending
=========

New base force fields
---------------------

Add them to ``_PRESETS`` in :mod:`unleap.amber`, as ``name -> (dat file, tleap preamble)``.
:func:`~unleap.amber.preset_bases` lists only the presets whose ``.dat`` exists in ``AMBERHOME``,
and ``--base auto`` scores them all.

Multi-residue or covalent ligands
---------------------------------

This would need:

* :func:`unleap.extract.strip` to accept several residues, and to stop rejecting bonds to
  outside atoms.
* :func:`unleap.lib.write` to emit connected units (head/tail atoms).
* A policy for bonded terms that cross the cut.

NBFIX
-----

ParmEd can write ``LJEDIT`` sections, but ``AmberParm`` doesn't populate ``atom_type.nbfix`` on
load. You'd have to build the pairs from the A/B coefficient tables found in
``extract._nonbonded_checks``.

Conventions
-----------

* User-facing failures raise :class:`unleap.UnleapError`. The CLI turns them into exit code 2
  with a one-line message.
* Warnings that still allow output go in ``ligand.warnings`` and appear as ``WARNING:`` lines in
  ``report.txt``.
* Any fix for a real topology needs a regression test that reproduces it synthetically.
