Output files
============

.. code-block:: text

   out/
     LIG.frcmod  LIG.full.frcmod  LIG.lib  LIG.mol2  LIG.sdf
     tleap.in                 # lines to load the outputs
     report.txt               # everything printed, incl. PASS/FAIL details and warnings
     LIG.orig.parm7/.rst7     # the ligand stripped from the input (non-periodic); what validation compares against
     check/                   # tleap/sander inputs, logs and rebuilt topologies from validation

``LIG`` is the ligand's residue name, or ``--name``.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - File
     - Contents
   * - ``LIG.frcmod``
     - Only the terms missing from, or different in, the base force field, like parmchk2 output.
       **Every improper is always listed**, as parmchk2 does (see :doc:`../dev/amber_notes`).
       Load it after the base, e.g. ``source leaprc.gaff2``.
   * - ``LIG.full.frcmod``
     - Every MASS / BOND / ANGLE / DIHE / IMPROPER / NONB term the ligand uses. Self-contained:
       rebuilds the ligand with no leaprc loaded.
   * - ``LIG.lib``
     - OFF library unit: atom names, types, exact charges, elements, connectivity and coordinates.
   * - ``LIG.mol2``
     - Tripos mol2 with GAFF atom types in the type column, exact charges (``USER_CHARGES``) and
       kekulized bond orders inferred by RDKit.
   * - ``LIG.sdf``
     - The same molecule with bond orders and formal charges, for RDKit and viewers. RDKit's mol2
       reader only accepts SYBYL atom types, so it cannot read any GAFF-typed mol2 (including
       antechamber's). Only written if bond orders were assigned.

Precision
---------

* Charges: 6 decimals, as in antechamber output. The prmtop round trip is exact.
* Bond ``req``: 7 decimals. Angle ``theteq``: 6 decimals, corrected for leap's degree conversion,
  so values come back as typed (113.02, not 113.0200487).
* Coordinates in the lib/mol2: about 1e-4 Å.
