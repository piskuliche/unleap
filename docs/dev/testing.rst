Testing
=======

Unit and integration tests
--------------------------

.. code-block:: bash

   python -m pytest -q          # ~5 s; skips everything if AmberTools isn't available

The fixtures are built per session in pytest's tmp dir (``tests/conftest.py``): acetaminophen
(GAFF2 and GAFF) and 4-chlorobenzoate next to an ff14SB peptide. Charges are antechamber
``-c gas`` for speed.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Area
     - Tests
   * - Round trip and diff ⊆ parmchk2
     - ``test_gaff2_roundtrip``
   * - Custom dihedral/LJ edits recovered exactly
     - ``test_custom_parameters_are_recovered``
   * - Auto-detect in a complex, charged-ligand bond orders
     - ``test_complex_autodetect_and_bond_orders``, ``test_sdf_is_rdkit_readable``
   * - HMR undo / ``--keep-masses`` refusal
     - ``test_hmr_is_undone``, ``test_keep_masses_rejects_split_types``
   * - Base auto-detect (gaff1)
     - ``test_gaff1_base_autodetect``
   * - NetCDF trajectory and REMD-style rst7 input
     - ``test_netcdf_trajectory_and_remd_rst7``
   * - Reversed improper order, leap angle correction
     - ``test_reversed_improper_order``, ``test_leap_angle_correction_recovers_input``
   * - Error paths
     - ``test_ambiguous_autodetect``, ``test_multi_residue_mask_rejected``, ``test_nbfix_warns``
   * - Validation is not vacuous
     - ``test_validation_catches_missing_terms`` (drops the dihedrals from the diff and expects FAIL)

When you fix a new real-world case, add a fixture that reproduces it. Items 3–5 in
:doc:`amber_notes` were invisible on the synthetic systems until real RBFE/ABFE topologies
exposed them.

Real-data sweep
---------------

.. code-block:: bash

   python scripts/vcm_sweep.py path/to/ABFE /tmp/vcm_sweep   # directory holding vcm_* jobs; ~10 min

``scripts/vcm_sweep.py`` runs unleap over every
``ABFE/vcm_*/results/.../t1_{complex,binder}_abfeproduction`` directory (its parm7 plus the λ=0
production rst7). It then rebuilds each job's ``inputs/abfe/<L>/binder_<L>.*`` with tleap and
compares that to the extracted ligand with :func:`~unleap.validate.compare_terms` and
:func:`~unleap.validate.compare_energies`.

That is the strongest end-to-end check available: it asks "is what we recovered what was put in?",
not just "do our files reproduce the topology?". It writes one directory per job and leg, plus
``summary.json``. It's worth re-running after any change to ``frcmod.py`` or ``validate.py``.

Last result: 52/52 topologies pass, and 46/46 with inputs are identical to their inputs.

Docs build
----------

.. code-block:: bash

   make -C docs html SPHINXOPTS=-W     # warnings are errors
