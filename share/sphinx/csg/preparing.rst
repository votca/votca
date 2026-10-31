.. _preparing:

Preparing coarse-grained runs
=============================

Preliminary note
~~~~~~~~~~~~~~~~

The coarse-grained run requires the molecule topology on the one hand
and suitable potentials on the other. In this chapter, the generation of
coarse-grained runs is described next, followed by a post-processing of
the potential.

If the potential is of such a form that it allows direct fitting of a
functional form, the section on post-processing can be skipped. Instead,
a program of choice should be used to fit a functional form to the
potential. Nevertheless, special attention should be paid to units
(angles, bond lengths). The resulting curve can then be specified in the
MD package used for simulation. However, most potentials don’t allow an
easy processing of this kind and tabulated potentials have to be used.

Generating a topology file for a coarse-grained run
---------------------------------------------------

**WARNING: This section describes experimental features. The exact names and options of the program might change in the near future. The section is specific to GROMACS support though a generalization for other MD packages is planned.**

The mapping definition is close to a topology needed for a coarse
grained run. To avoid redundant work, ``csg_gmxtopol`` can be used to automatically
generate a gromacs topology based on an atomistic reference system and a
mapping file.

At the current state, ``csg_gmxtopol`` can only generate the topology for the first
molecule in the system. If more molecule types are present, a special
tpr file has to be prepared. The program can be executed by

.. code:: bash

      csg_gmxtopol --top topol.tpr --cg map.xml --out cgtop

which will create a file ``cgtop.top``. This file includes the topology
for the first molecule including definitions for atoms, bonds, angles
and dihedrals. It can directly be used as a topology in GROMACS, however
the force field definitions (atom types, bond types, etc.) still have to
be added manually.

.. _preparing_post-processing_of_the_potential:

Post-processing of the potential
--------------------------------

VOTCA no longer provides tools for post-processing potentials (the former
``csg_resample`` and ``csg_call`` table scripts were removed); these steps
have to be done with other tools. They roughly consist of

-  clipping poorly sampled (border) regions, e.g. by deleting them from the
   table or by setting their flag to ``o`` (out of range),

-  resampling the potential to the grid required by the MD package,

-  extrapolating the potential into the clipped regions,

-  exporting the table in the format of the MD package.

Tables written by ``csg_stat`` and ``csg_fmatch`` have the columns
``x y [error] flag``, where the flag ``i`` marks values inside the sampled
range and ``o`` values outside of it.
