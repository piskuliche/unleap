Validation
==========

Unless you pass ``--no-check``, unleap rebuilds the ligand with tleap three ways and compares each
rebuild with the stripped original, ``LIG.orig.parm7``:

1. ``lib + diff frcmod``, with the base leaprc
2. ``lib + full frcmod``, with no leaprc
3. ``mol2 + diff frcmod``

Each comparison has two parts:

**Term by term**, keyed by atom indices:
bonds, angles, multi-term dihedrals (including SCEE/SCNB), impropers, charges, LJ, masses
(after the HMR undo), atom names and types. A difference is reported with atom names, e.g.
``improper O6-C6-N1-C5: missing in rebuild``.

**By energy**:
a sander vacuum single point (``ntb=0``, ``cut=999``) of BOND, ANGLE, DIHED, VDWAALS, EEL,
1-4 VDW and 1-4 EEL, which must agree to within 2e-3 kcal/mol. The rebuilt topology is scored
on the *original* coordinates, so only parameters are being compared.

If tleap itself fails, for example on a missing parameter, that also counts as a FAIL, and
``report.txt`` gives the path to the tleap log.

Reading a failure
-----------------

.. code-block:: text

   round-trip validation (vs stripped original):
     [FAIL] lib + diff frcmod
            improper O6-C6-C5-N1: extra in rebuild
            improper O6-C6-N1-C5: missing in rebuild
     [PASS] lib + full frcmod

* If **only the diff rebuilds fail**, the base force field reproduces something differently from the
  original build. Use ``LIG.full.frcmod``, or try a different ``--base``.
* If **all rebuilds fail**, unleap doesn't handle something in this topology. Please report it with
  the inputs.

Tested on
---------

* Synthetic GAFF/GAFF2 systems (water box, peptide complex, custom dihedral/LJ edits, HMR).
* Production free-energy topologies from an RNA-ligand campaign: RBFE dual topologies
  (``unisc.parm7``, both end-state ligands), ABFE trajectories, and 26 amberstudio ABFE jobs
  (complex and binder legs: 52 HMR topologies). All pass. For every job that kept its ``inputs/``,
  the recovered parameters are identical to the job's original ``binder_*.mol2/.frcmod/.lib``.
