Diff frcmod matching rules
==========================

:func:`unleap.frcmod.compare_to_base` decides, term by term, whether the base force field already
reproduces the ligand's parameter. It mirrors tleap's lookup order:

.. list-table::
   :header-rows: 1
   :widths: 14 40 46

   * - Term
     - Lookup
     - Equal if
   * - atom type
     - by name
     - mass, rmin and ε within ``_TOL`` (rel/abs 1e-5)
   * - bond / angle
     - explicit, either direction
     - k and eq within ``_TOL``
   * - dihedral
     - explicit (either direction), then ``X-b-c-X``
     - same term list after dropping zero-k terms: periodicity, phase (mod 360), k, SCEE, SCNB
   * - improper
     - central third; outer atoms permuted; 0, then 1, then 2 ``X`` wildcards
     - k, periodicity, phase (no 1-4 scaling)

A term that is missing from the base or differs is written explicitly. tleap gives explicit
frcmod entries loaded after the base priority over the base, which reproduces the original
build. Impropers are all written regardless (see :doc:`amber_notes`, item 5).

Base scoring
------------

:func:`unleap.frcmod.score_base` is the fraction of the ligand's bonded terms plus LJ types that
a base reproduces. Masses are ignored, because HMR may have altered them. ``--base auto`` takes
the best preset from :func:`unleap.amber.preset_bases`.

What the diff cannot express
----------------------------

A frcmod can add or override terms, but it cannot remove one. If the base would add a term the
original build didn't have, validation reports ``extra in rebuild``. ``LIG.full.frcmod`` does not
have this problem.
