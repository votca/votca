Force matching
==============

.. figure:: usage/fig/flow_fmatch.eps
   :alt: Flowchart to perform force matching.
   :width: 50.0%

   Flowchart to perform force matching.

The force matching algorithm with cubic spline basis is implemented in
the utility. A list of available options can be found in the reference
section of (command ``–h``).

Program input
-------------

needs an atomistic reference run to perform coarse-graining. Therefore,
the trajectory file *must contain forces* (note that there is a suitable
option in the ``.mdp`` file), otherwise will not be able to run.

In addition, a mapping scheme has to be created, which defines the
coarse-grained model (see ). At last, a control file has to be created,
which contains all the information for coarse-graining the interactions
and parameters for the force-matching run. This file is specified by the
tag in the format. An example might look like the following

::

       <cg>
         <!--fmatch section -->
         <fmatch>
           <!--Number of frames for block averaging -->
           <frames_per_block>6</frames_per_block>
           <!--Constrained least squares?-->
           <constrainedLS>false</constrainedLS>
         </fmatch>
         <!-- example for a non-bonded interaction entry -->
         <non-bonded>
           <!-- name of the interaction -->
           <name>CG-CG</name>
           <type1>A</type1>
           <type2>A</type2>
           <!-- fmatch specific stuff -->
           <fmatch>
             <min>0.27</min>
             <max>1.2</max>
             <step>0.02</step>
             <out_step>0.005</out_step>
           </fmatch>
         </non-bonded>
       </cg>

Similarly to the case of spline fitting (see sec. [sec:ref\_programs] on
), the parameters ``min`` and ``max`` have to be chosen in such a way as
to avoid empty bins within the grid. Determining ``min`` and ``max`` by
using is recommended (see sec. [sec:setting\_files]). A full description
of all available options can be found in sec. [sec:ref\_options].

Program output
--------------

produces a separate ``.force`` file for each interaction, specified in
the CG-options file (option ). These files have 4 columns containing
distance, corresponding force, a table flag and the force error, which
is estimated via a block-averaging procedure. If you are working with an
angle, then the first column will contain the corresponding angle in
radians.

To get table-files for , integrate the forces in order to get potentials
and do extrapolation and potentially smoothing afterwards.

Output files are not only produced at the end of the program execution,
but also after every successful processing of each block. The user is
free to have a look at the output files and decide to stop , provided
the force error is small enough.

Integration and extrapolation of .force files 
----------------------------------------------

To convert forces (``.force``) to potentials (``.pot``), tables have to
be integrated. To use the built-in integration command from the
scripting framework, execute

::

     $csg_call table integrate CG-CG.force minus_CG-CG.pot
     $csg_call table linearop minus_CG-CG.d CG-CG.d -1 0

This command calls the script, which integrates the force and writes the
potential to the ``.pot`` file.

In general, each potential contains regions which are not sampled. In
this case or in the case of further post-processing, the potential can
be refined by employing resampling or extrapolating methods. See sec.
[sec:post\_processing] for further details.
