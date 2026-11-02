unleap
======

**Get a ligand's parameter files back out of a built Amber system.**

tleap takes a ligand's ``mol2`` / ``frcmod`` / ``lib`` and builds a system (``parm7`` + ``rst7``).
``unleap`` goes the other way. Give it a complex or solvated system and it strips the ligand out,
regenerates those files, and then proves they are right. It rebuilds the ligand with tleap and
checks every parameter and a sander energy breakdown against the original topology.

.. code-block:: bash

   unleap -p complex.parm7 -c complex.rst7 -o lig_out

or, from Python, ``unleap.extract_ligand("complex.parm7", "complex.rst7")`` (:doc:`user/python_api`).

Typical uses:

* The original ligand inputs are gone, or you no longer trust they match what was simulated.
* You want the exact parameters used in a production run (RBFE/ABFE/VCM job) to reuse elsewhere.
* You want to check that a run's topology really used the ligand files you think it did.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   user/installation
   user/quickstart
   user/options
   user/python_api
   user/outputs
   user/validation
   user/limitations

.. toctree::
   :maxdepth: 2
   :caption: Developer guide

   dev/architecture
   dev/amber_notes
   dev/diff_rules
   dev/testing
   dev/extending

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   changelog
