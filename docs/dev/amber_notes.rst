Amber and ParmEd behaviors
==========================

Each of these was found on real data and has a regression test (see :doc:`testing`).
Don't "simplify" these away.

1. Structure slicing breaks AmberParm
-------------------------------------

``struct[mask]`` keeps the full LJ tables, so the saved parm7 is malformed::

   FLAG LENNARD_JONES_ACOEF has 66 elements; expected 45

``Structure.copy(AmberParm)`` drops ``parm_data`` entirely. Reload the file and call
``AmberParm.strip()``, which rebuilds the tables (:func:`unleap.extract.strip`).

2. ParmEd's frcmod writer loses precision
-----------------------------------------

``AmberParameterSet.write`` formats ``req``, ``theteq`` and masses with 3 decimals, but GAFF 2.2
uses 4-decimal bond lengths (1.3986). It also re-sorts improper atoms alphabetically.
:func:`unleap.frcmod.write` is our own writer in the same layout: 7 decimals for ``req``,
6 for ``theteq``, and impropers kept in topology order.

3. leap's DEGTORAD is 0.0174533
-------------------------------

``leap/src/leap/basics.h`` defines ``DEGTORAD 0.0174533``, not π/180. Every angle and phase tleap
writes is inflated by about 4.3e-7 when read back exactly, so 113.02 becomes 113.0200487.
Writing that back inflates it again. :func:`unleap.frcmod.full_params` divides by
:data:`unleap.frcmod.LEAP_DEG_CORRECTION`, but only when the corrected values are "rounder"
(``_leap_inflated``), so a topology written by another tool is not over-corrected.

4. tleap may store an improper reversed
---------------------------------------

The prmtop atom order can be ``l-k-j-i``, which puts the central atom second. It is the same
torsion angle, but "central = third atom" is then wrong. For one production ligand this turned the original
``cc-ns-c -o`` into ``o -c -ns-cc``, which tleap then never matched.
:func:`unleap.frcmod.improper_atoms` finds the central atom from bonding (bonded to the other
three) and reorders. The validation keys use the same function.

5. Wildcard and explicit impropers give different atom orders
-------------------------------------------------------------

If tleap matches ``X -X -c -o`` rather than an explicit ``cc-ns-c -o``, it stores the four atoms
in a different order. The improper angle, and so the energy, changes slightly even though the
parameters are identical. So the diff frcmod **always lists every improper explicitly**, as
parmchk2 does. :func:`unleap.frcmod.lookup_improper` only feeds the base score.

6. from_structure keys impropers under every permutation
--------------------------------------------------------

``AmberParameterSet.from_structure`` stores each improper under every outer-atom permutation.
unleap dedupes on ``(central type, sorted outer types)``.

7. tleap requires MASS and NONBON together
------------------------------------------

A frcmod with NONBON but no MASS fails::

   Modified force field files must have both MASS and NONBOND entries or neither.

The writer always emits both for every atom type in the set.

8. Coordinate reading
---------------------

ParmEd's format sniffing rejects ASCII restarts whose natom line also holds a temperature
(pmemd/REMD: ``52598  0.1400000E+05  0.2980000E+03``). :func:`unleap.extract.read_coordinates`
uses ``Rst7.open`` directly. For NetCDF it uses ``NetCDFTraj`` (a ``frame`` dimension exists) or
``NetCDFRestart``, and reads a single frame from ``_ncfile.variables``, so it never loads a 750 MB
trajectory just to take one frame.

9. The lib stores coordinates at about 1e-4 Å
---------------------------------------------

For stiff angles that alone moves the energy by about 1e-3 kcal/mol, so validation scores each
rebuilt topology on the *original* coordinates. The atom order is identical.

10. RDKit can't read GAFF-typed mol2
------------------------------------

RDKit's mol2 reader expects SYBYL types (``Element 'o' not found``). That's why the SDF is
written too.
