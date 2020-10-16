Preparing coarse-grained runs
=============================

Preliminary note
~~~~~~~~~~~~~~~~

The coarse-grained run requires the molecule topology on the one hand
and suitable potentials on the other. In this chapter, the generation of
coarse-grained runs is decribed next, followed by a post-processing of
the potential.

If the potential is of such a form that it allows direct fitting of a
functional form, the section on post-processing can be skipped. Instead,
a program of choice should be used to fit a functional form to the
potential. Nevertheless, special attention should be paid to units
(angles, bondlengths). The resulting curve can then be specified in the
MD package used for simulation. However, most potentials don’t allow an
easy processing of this kind and tabulated potentials have to be used.

Generating a topology file for a coarse-grained run
---------------------------------------------------

WARNING: This section describes experimental features. The exact names and options of the program might change in the near future. The section is specific to support though a generalization for other MD packages is planned.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mapping definition is close to a topology needed for a coarse
grained run. To avoid redundant work, can be used to automatically
generate a gromacs topology based on an atomistic reference system and a
mapping file.

At the current state, can only generate the topology for the first
molecule in the system. If more molecule types are present, a special
tpr file has to be prepared. The program can be executed by

::

      csg_gmxtopol --top topol.tpr --cg map.xml --out cgtop

which will create a file ``cgtop.top``. This file includes the topology
for the first molecule including definitions for atoms, bonds, angles
and dihedrals. It can directly be used as a topology in , however the
force field definitions (atom types, bond types, etc.) still have to be
added manually.

Post-processing of the potential
--------------------------------

The package provides a collection of scripts to handle potentials. They
can be modified, refined, integrated or inter- and extrapolated. These
scripts are the same ones as those used for iterative methods in chapter
[sec:iterative\_methods]. Scripts are called by . A complete list of
available scripts can be found in sec. [sec:csg\_table].

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
first. This is simply done by editing the file and by deleting those
values.

Alternatively, manually check the range where the potential still looks
good and is not to noisy and set the flags in the potential file of the
bad parts by hand to ``o`` (for ``out of range``). Those values will
later be extrapolated and overwritten.

Resampling
~~~~~~~~~~

Use the command

::

      csg_resample --in table.pot --out table_resample.pot \
                   --grid min:step:max

to resample the potential given in file – from to with a grid spacing of
steps. The result is written to the file specified by . Additionally,
allows the specification of spline interpolation (), the calculation of
derivatives () and comments (). Check the help () for further
information.

It is important to note that the values and *don’t* correspond to the
minimum and maximum value in the input file, but to the range of values
the potential is desired to cover after extrapolation. Therefore, values
in :math:`[ \min,\max ]` that are not covered in the file are
automatically marked by a flag ``o`` (for ``out of range``) for
extrapolation in the next step.

The potential *don’t* have to start at 0, this is done by the export
script (to xvg) automatically.

Extrapolation
~~~~~~~~~~~~~

The following line

::

      csg_call table extrapolate [options] table_resample.pot \
               table_extrapolate.pot

calls the extrapolation procedure, which processes the range of values
marked by . The input file is created in the last step.

After resampling, all values in the potential file that should be used
as a basis for extrapolation are marked with an ``i``, while all values
that need extrapolation are marked by ``o``. The command above now
extrapolates all ``o`` values from the ``i`` values in the file.
Available options include averaging over a certain number of points (),
changing the functional form (, default is quadratic), extrapolating
just the left or right region of the file () and setting the curvature
().

The output of the extrapolation step can now be used for the
coarse-grained run. If is used as a molecule dynamics package, the
potential has to be converted and exported to a suitable format as
described in the final step.

Exporting the table
~~~~~~~~~~~~~~~~~~~

Finally, the table is exported to . The conversion procedure requires a
small xml file as shown below:

::

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

where is the and is just a number slightly smaller than the upper value
of single/ double precision. The value given in corresponds to the value
of .

Using the file above, call

::

      csg_call --options table.xml --ia-type non-bonded --ia-name XXX \
        convert_potential gromacs table_extrapolate.pot table.xvg

to convert the extrapolated values in to (The file will contain the C12
parts only which are stored in the sixth und seventh column, this can be
changed by adding the option (for the fourth and fiveth column) or
option (for the second and third column) after . Ensure compatibility
with the topology. See the manual for further information).

To obtain a bond table, run

::

      csg_call --ia-type bond --ia-name XXX --options table.xml \
      convert_potential gromacs table_extrapolate.pot table.xvg

It is also possible to use ``angle`` and ``dihedral`` as type as well,
but make to sure to have a ``bonded`` section similar to the
``non-bonded`` section above with the corresponding interaction name.

Internally will do the following steps:

-  Resampling of the potential from 0 (or -180 for dihedrals) to
   ``table_end`` (or 180 for angles and dihedrals) with step size
   ``table_bins``. This is needed for gromacs the table must start with
   0 or -180.

-  Extrapolate the left side (to 0 or -180) exponentially

-  Extrapolate the right side (to ``table_end`` or 180) exponentially
   (or constant for non-bonded interactions)

-  Shift it so that the potential is zero at ``table_end`` for
   non-bonded interactions or zero at the minimum for bonded interaction

-  Calculate the force (assume periodicity for dihedral potentials)

-  Write to the format needed by gromacs

An example on non-bonded interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

      csg_call pot shift_nonbonded table.pot table.pot.refined
      csg_resample --grid 0.3:0.05:2 --in table.pot.refined \
               --out table.pot.refined
      csg_call table extrapolate --function quadratic --region left \
               table.pot.refined table.pot.refined
      csg_call table extrapolate --function constant --region right \
               table.pot.refined table.pot.refined

Alternatives
------------

Additionally to the two methods described above, namely (a) providing
the MD package directly with a functional form fitted with a program of
choice or (b) using , and , another method would be suitable. This is
integrating the force table as follows

::

      -Integrate the table
      $csg_call table integrate force.d minus_pot.d
      -multiply by -1
      $csg_call table linearop minus_pot.d pot.d -1 0
