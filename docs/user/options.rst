Options
=======

.. code-block:: text

   unleap -p PARM -c CRD [--frame N] [--mask MASK | --resname NAME] [-o OUTDIR]
          [--name NAME] [--base BASE] [--keep-masses] [--no-check]

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Option
     - Meaning
   * - ``-p, --parm``
     - Input topology (``parm7``/``prmtop``). Required.
   * - ``-c, --crd``
     - Coordinates: ASCII ``rst7``/``inpcrd`` (including pmemd/REMD restarts with a temperature
       in the header), NetCDF restart, or NetCDF trajectory. Required.
   * - ``--frame N``
     - Frame to use from a NetCDF trajectory. 0-based; ``-1`` = last. Default ``0``.
   * - ``--mask MASK``
     - Amber mask selecting the ligand residue, e.g. ``':LIG'`` or ``':5'``.
   * - ``--resname NAME``
     - Select the ligand by residue name (shorthand for ``--mask :NAME``).
   * - ``-o, --outdir``
     - Output directory. Default ``unleap_out``.
   * - ``--name NAME``
     - Rename the residue/unit in the outputs. Default: keep the original name.
   * - ``--base BASE``
     - Base force field for the diff frcmod: ``auto`` (default), ``gaff``, ``gaff2``,
       ``gaff211``, or a path to a ``.dat`` file.
   * - ``--keep-masses``
     - Do not undo hydrogen mass repartitioning.
   * - ``--no-check``
     - Skip the tleap/sander round-trip validation.

Selecting the ligand
--------------------

With neither ``--mask`` nor ``--resname``, unleap picks the single residue that is not a standard
amino acid (including protonation variants and caps), nucleotide, water or ion. If there is more
than one candidate, for example the two ligands in an RBFE dual topology, it stops and lists them.

The selection must be **exactly one whole residue** that is **not covalently bonded** to anything
outside the selection.

Choosing the base force field
-----------------------------

``--base auto`` scores gaff, gaff2 (2.2) and gaff211 by how many of the ligand's terms each one
reproduces, and uses the best. The first line of the diff frcmod records which base it was
written against and the tleap line that loads it.

.. list-table::
   :header-rows: 1

   * - Base
     - tleap line
   * - ``gaff``
     - ``source leaprc.gaff``
   * - ``gaff2``
     - ``source leaprc.gaff2``
   * - ``gaff211``, or a ``.dat`` path
     - ``loadamberparams <file>.dat``

Pass ``--base`` explicitly if, for example, the system was built against an older gaff2 that you
have as a ``.dat`` file.

Hydrogen mass repartitioning
----------------------------

Production topologies often have HMR applied (H masses ≈ 3.024). unleap detects this and restores
the original masses. Each H gets its type's standard mass, and the excess goes back to the heavy
atom it is bonded to, which reverses HMR exactly.

With ``--keep-masses`` the topology's masses are kept as they are. That usually fails on HMR
systems: a ``ca`` carbon bonded to an H and one without an H end up with different masses, and
one frcmod atom type cannot hold both.

Exit codes
----------

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Code
     - Meaning
   * - ``0``
     - Files written; all validation checks passed (or ``--no-check``).
   * - ``1``
     - Files written, but at least one validation check failed. Details in ``report.txt``.
   * - ``2``
     - Input or setup error (message on stderr); nothing reliable was written.
