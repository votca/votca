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
or can be output directly by calling

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
the sampling package. The script runs from the package only on one cpu.
Our task will be to change the script so that uses 8 cpus, which is
basically the same as adding mpirun options in .

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
clear, one can try to find the right script in . Alternatively, check
the folder

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
in the setting file (see ).

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
