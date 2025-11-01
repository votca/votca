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

The post-processing roughly consists of the following steps (see further
explanations below):

-  (manually) clipping poorly sampled (border) regions

-  resampling the potential in order to change the grid to the proper
   format ()

-  extrapolation of the potential at the borders ( table extrapolate)

-  exporting the table to xvg ( convert\_potential gromacs)

Clipping of poorly sampled regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Regions with an irregular distribution of samples should be deleted
first. This is simply done by editing the ``.pot`` file and by deleting
those values.

Alternatively, manually check the range where the potential still looks
good and is not to noisy and set the flags in the potential file of the
bad parts by hand to ``o`` (for ``out of range``). Those values will
later be extrapolated and overwritten.

Resampling
~~~~~~~~~~

Use the command

.. code:: bash

      csg_resample --in table.pot --out table_resample.pot \
                   --grid min:step:max

to resample the potential given in file –``table.pot`` from ``min`` to
``max`` with a grid spacing of ``step`` steps. The result is written to
the file specified by ``out``. Additionally, ``csg_resample`` allows the specification of
spline interpolation (``spfit``), the calculation of derivatives
(``derivative``) and comments (``comment``). Check the help (``help``)
for further information.

It is important to note that the values ``min`` and ``max`` *don’t*
correspond to the minimum and maximum value in the input file, but to
the range of values the potential is desired to cover after
extrapolation. Therefore, values in :math:`[ \min,\max ]` that are not
covered in the file are automatically marked by a flag ``o`` (for
``out of range``) for extrapolation in the next step.

The potential *don’t* have to start at 0, this is done by the export
script (to xvg) automatically.

Extrapolation
~~~~~~~~~~~~~

After resampling, all values in the potential file that should be used
as a basis for extrapolation are marked with an ``i``, while all values
that need extrapolation are marked by ``o``. The command above now
extrapolates all ``o`` values from the ``i`` values in the file.
Available options include averaging over a certain number of points
(``avgpoints``), changing the functional form (``function``, default is
quadratic), extrapolating just the left or right region of the file
(``region``) and setting the curvature (``curvature``).

The output ``table_extrapolate.pot`` of the extrapolation step can now
be used for the coarse-grained run. If GROMACS is used as a molecule
dynamics package, the potential has to be converted and exported to a
suitable GROMACS format as described in the final step.

Exporting the table
~~~~~~~~~~~~~~~~~~~

Finally, the table is exported to ``xvg``. The conversion procedure
requires a small xml file ``table.xml`` as shown below:

.. code:: xml

      <cg>
        <non-bonded>
          <name>XXX</name>
          <step>0.01</step>
        </non-bonded>
        <inverse>
          <gromacs>
            <pot_max>1e8</pot_max>
            <table_end>8.0</table_end>
            <table_bins>0.002</table_bins>
          </gromacs>
        </inverse>
      </cg>

where ``<table_end>`` is the GROMACS ``rvdw+table_extension`` and
``<pot_max>`` is just a number slightly smaller than the upper value of
single/ double precision. The value given in ``<table_bins>``
corresponds to the ``step`` value of
``csg_resample -grid min:step:max``.


