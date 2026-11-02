Limitations and troubleshooting
===============================

Limitations
-----------

* **One residue only, non-covalent only.** Multi-residue or covalently attached ligands are refused.
* **NBFIX / changeLJPair edits and 12-6-4 C coefficients** on the ligand's atom types cannot be
  written to a frcmod. unleap warns, and they are left out.
* **CMAP terms** on ligand atoms are dropped, with a warning.
* **One atom type name, one set of parameters.** If two atoms share a type name but have different
  LJ parameters or masses, unleap stops, because one frcmod entry cannot describe both.
* **Bond orders are inferred** by RDKit from connectivity, coordinates and the net charge. If that
  fails, the mol2 gets all single bonds (tleap does not care), no SDF is written, and you get a
  warning.
* **Amber force fields only.** CHARMM (chamber) and AMOEBA topologies are rejected.
* **The diff frcmod cannot remove a term** that the base force field would add. Validation reports
  such a case as ``extra in rebuild``. Use the full frcmod.

Troubleshooting
---------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Message
     - Meaning / fix
   * - ``AMBERHOME is not set``
     - Source AmberTools first (``source $AMBERHOME/amber.sh``, or activate a conda env with ``ambertools``).
   * - ``auto-detect is ambiguous; candidates: ...``
     - Pass ``--resname`` or ``--mask``.
   * - ``mask ... spans N residues``
     - Select a single residue.
   * - ``... bonded outside the selection (covalent ligand)``
     - Not supported.
   * - ``atom type 'X' has inconsistent mass/LJ ... (hydrogen mass repartitioning?)``
     - Drop ``--keep-masses``.
   * - ``--frame N out of range``
     - The trajectory has fewer frames. Negative values count from the end.
   * - ``could not read coordinates from ...``
     - The file isn't an rst7/inpcrd or NetCDF file, or its atom count differs from the topology.
   * - ``[FAIL] ...`` lines
     - See :doc:`validation`. Look at ``report.txt`` and ``check/*.log``.
