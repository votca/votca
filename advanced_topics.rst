Advanced topics
===============

Customization
-------------

Each sub-step of an iteration and all direct calls can be adjusted to
the user needs. The internal part of the iterative framework is
organized as follows: all scripts are called using two keywords

::

      csg_call key1 key2

For example, ``csg_call update imc`` calls the ``update`` script for the
inverse Monte Carlo procedure. The corresponding keywords are listed in
sec. [sec:csg\_table] or can be output directly by calling

::

      csg_call --list

It is advised not to change already implemented scripts. To customize a
script or add a new one, copy the script to your own directory (set by )
and redirect its call by creating your own ``csg_table`` file in this
directory which looks like this

::

      key1 key2 script1 options
      key3 key4 script2

If the local keys are already in use, the existing call will be
overloaded.

As an example, we will illustrate how to overload the script which calls
the sampling package. The script runs from the GROMACSpackage only on
one cpu. Our task will be to change the script so that GROMACSuses 8
cpus, which is basically the same as adding mpirun options in .

First we find out which script calls :

::

      csg_call --list | grep gromacs

The output should look as follows

::

      init gromacs initalize_gromacs.sh
      prepare gromacs prepare_gromacs.sh
      run gromacs run_gromacs.sh
      pressure gromacs calc_pressure_gromacs.sh
      rdf gromacs calc_rdf_gromacs.sh
      imc_stat gromacs imc_stat_generic.sh
      convert_potential gromacs potential_to_gromacs.sh

the third line indicates the script we need. If the output of is not
clear, one can try to find the right script in sec. [sec:csg\_table].
Alternatively, check the folder

::

      <csg-installation>/share/scripts/inverse

for all available scripts.

Analyzing the output of

::

      csg_call --cat run gromacs

we can conclude that this is indeed the script we need as the content
(in shorted form is):

::

      critical mdrun

Now we can create our own ``SCRIPTDIR``, add a new script there, make it
executable and overload the call of the script:

::

      mkdir -p SCRIPTDIR
      cp `csg_call --quiet --show run gromacs` SCRIPTDIR/my_run_gromacs.sh
      chmod 755 SCRIPTDIR/my_run_gromacs.sh
      echo "run gromacs my_run_gromacs.sh" >> SCRIPTDIR/csg_table

Please note that ``my_run_gromacs.sh`` is the name of the script and
``SCRIPTDIR`` is the custom script directory, which can be a global or a
local path. Now we change the last line of ``my_run_gromacs.sh`` to:

::

      critical mpirun -np 8 mdrun

This completes the customization. Do not forget to add ``SCRIPTDIR`` to
in the setting XMLfile (see sec. [sec:ref\_options]).

You can check the new script by running:

::

      csg_call --scriptdir SCRIPTDIR --list
      csg_call --scriptdir SCRIPTDIR --run run gromacs

Finally, do not forget to remove the license infomation and change the
version number of the script.

Used external packages
----------------------

GroMaCS
~~~~~~~

Get it from

-  
-  

ESPResSo
~~~~~~~~

Get it from

DL\_POLY
~~~~~~~~

Get it from

Gnuplot
~~~~~~~

Get it from

LAMMPS
~~~~~~

Get it from


DL\_POLY interface
------------------

**WARNING: The DL\_POLYinterface is still experimental (in development)
but it does support the Iterative Boltzmann Inversion and Inverse Monte
Carlo schemes. The Force Matching might work as well, although it has
not been tested thoroughly.**


The DL\_POLYinterface fully supports coarse-grain mapping of a full-atom
system previuosly simulated with any version of DL\_POLY, including
DL\_POLY-Classic. However, the full optimization of the effective
potentials with the aid of iterative methods will only become possible
when the new release of DL\_POLY-4(4.06) is made public; the reason
being the incapability of earlier DL\_POLYversions of using
user-specified tabulated force-fields for intramolecular, aka “bonded”,
interactions: bonds, angles, dihedral angles (torsions). Below the
coarse-graining and CG force-field optimization with the aid of the
latest DL\_POLY-4version (4.06+) are outlined.

Running votcawith DL\_POLY-4as MD simulation engine is very similar to
doing so with GROMACS. The three types of required input files in the
case of DL\_POLYare: CONTROL – containing the simulation directives and
parameters (instead of ``.mdp`` file for GROMACS), FIELD – the topology
and force-field specifications (instead of ``.top`` and ``.tpr`` files),
and CONFIG (instead of ``.gro`` file) – the initial configuration file,
containing the MD cell matrix and particle coordinates (it can also
include initial velocities and/or forces); for details see
DL\_POLY-4manual. Most of the votcatools and scripts described above in
the case of using GROMACSwill work in the same manner, with the
following conventional substitutions for the (default) file names used
in options for votcascripts, as necessary:

::

    .dlpf = the topology read from FIELD or written to FIELD_CGV
    .dlpc = the configuration read from CONFIG or written to CONFIG_CGV
    .dlph = the trajectory read from HISTORY or written to HISTORY_CGV

It is also possible to specify file names different from the standard
DL\_POLYconvention, in which case the user has to use the corresponding
dot-preceded extension(s); for example: FA-FIELD.dlpf instead of FIELD
or CG-HISTORY.dlph instead of HISTORY\_CGV (see section
[sec:ref\_programs], as well as the man pages or output of
votcacommands, with option ``—help``).

votcafollows the DL\_POLYconventions for file names and formats. Thus,
``csg_dlptopol`` and ``csg_map`` produce the CG topology (FIELD\_CGV by
default), configuration (CONFIG\_CGV), and/or trajectory (HISTORY\_CGV)
files fully compatible with and usable by DL\_POLY. **Note that the
ability of these tools to read and write a plethora of different file
formats provides means to convert input and output files between the
simulation packages supported by votca, e.g. GROMACS– DL\_POLYor vice
versa. The user is, however, strongly advised to check the resulting
files for consistency before using them).**

Similarly, the distribution analysis and potential/force generation
utilities, such as ``csg_stat`` and votcascripts, will read and write
DL\_POLY-formatted files; in particular, the tabulated force-field files
containing the potential and force/virial data: TABLE – for short-range
(VdW) “non-bonded” interactions, TABBND, TABANG and TABDIH – for
“bonded” interations: bonds, bending angles and dihedrals,
correspondingly (for the format details see DL\_POLY-4manual). Note,
however, that the latter three files can only be used by
DL\_POLY-4(4.06+).

The user is advised to search for “dlpoly” through the
``csg_defaults.xml``, ``csg_table`` files and in scripts located in
``share/votca/scripts/inverse/`` in order to find out about the xml-tags
and options specific for DL\_POLY; see also sections [sec:ref\_options]
and [sec:csg\_table].
