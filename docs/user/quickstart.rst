Quick start
===========

Single ligand, auto-detected
----------------------------

If the system has exactly one non-standard residue, nothing else is needed:

.. code-block:: bash

   unleap -p complex.parm7 -c complex.rst7 -o lig_out

RBFE dual topology
------------------

A ``unisc.parm7`` holds both end-state ligands, so pick one:

.. code-block:: bash

   unleap -p unisc.parm7 -c 0.00000000_ti.rst7 --resname L00 -o L00_out
   unleap -p unisc.parm7 -c 0.00000000_ti.rst7 --resname L04 -o L04_out

Straight from a trajectory
--------------------------

NetCDF trajectories work directly. ``--frame -1`` takes the last frame:

.. code-block:: bash

   unleap -p complex_LIG.parm7 -c complex_LIG_0.00000.nc --frame -1 -o LIG_out

amberstudio ABFE job
--------------------

Use the trial-1 production directory: its ``parm7`` plus the λ=0 production ``rst7``
(the fully interacting state):

.. code-block:: bash

   T=results/abfe_run/abfe/attempt-0001/LIG/t1_complex_abfeproduction
   unleap -p $T/complex_LIG.parm7 -c $T/complex_LIG_0.00000.rst7 -o LIG_out

Example output
--------------

.. code-block:: text

   unleap 0.1.0
   ligand: mask :L00 -> residue L00, 21 atoms, net charge -0.0000
   base force field: gaff2  (match: gaff 28.4%, gaff2 72.6%, gaff211 41.1%)
   HMR detected: H masses 3.024 -> 1.008, excess returned to the bonded heavy atoms
   diff frcmod: 0 atom, 3 bond, 8 angle, 10 dihedral, 8 improper
   wrote .../L00_out/L00.sdf
   wrote .../L00_out/L00.mol2
   ...
   round-trip validation (vs stripped original):
     [PASS] lib + diff frcmod
     [PASS] lib + full frcmod
     [PASS] mol2 + diff frcmod

Using the results
-----------------

.. code-block:: text

   source leaprc.gaff2
   loadamberparams L00.frcmod
   loadoff L00.lib

``tleap.in`` in the output directory contains exactly these lines. If you want one file that
works without any leaprc, load ``L00.full.frcmod`` instead.
