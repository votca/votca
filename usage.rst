Boltzmann Inversion
===================

provides a potential of mean force for a given degree of freedom.

t7cm |image|

It is mostly used for deriving *bonded* interactions from canonical
sampling of a single molecule in vacuum, e. g. for polymer
coarse-graining, where it is difficult to separate bonded and non-bonded
degrees of freedom :raw-latex:`\cite{Tschoep:1998}`. The non-bonded
potentials can then be obtained by using iterative methods or force
matching.

The main tool which can be used to calculate histograms, cross-correlate
coarse-grained variables, create exclusion lists, as well as prepare
tabulated potentials for coarse-grained simulations is . It parses the
whole trajectory and stores all information on bonded interactions in
memory, which is useful for interactive analysis. For big systems,
however, one can run out of memory. In this case can be used which,
however, has a limited number of tasks it can perform (see sec.
[sec:setting\_files] for an example on its usage).

Another useful tool is . It can be used to convert an atomistic
trajectory to a coarse-grained one, as it is discussed in
sec. [sec:trajectory].

To use one has to first define a mapping scheme. This is outlined in .
Once the mapping scheme is specified, it is possible to generate an
exclusion list for the proper sampling of the atomistic resolution
system.

Generating exclusion lists
--------------------------

Exclusion lists are useful when sampling from a special reference system
is needed, for example for polymer coarse-graining with a separation of
bonded and non-bonded degrees of freedom.

To generate an exclusion list, an atomistic topology without exclusions
and a mapping scheme have to be prepared first. Once the .tpr topology
and .xml mapping files are ready, simply run

::

      csg_boltzmann --top topol.tpr --cg mapping.xml --excl exclusions.txt

This will create a list of exclusions for all interactions that are not
within a bonded interaction of the coarse-grained sub-bead. As an
example, consider coarse-graining of a linear chain of three beads which
are only connected by bonds. In this case, will create exclusions for
all non-bonded interactions of atoms in the first bead with atoms of the
3rd bead as these would contribute only to the non-bonded interaction
potential. Note that will only create the exclusion list for the fist
molecule in the topology.

To add the exclusions to the topology of the molecule, either include
the file specified by the –excl option into the .top file as follows

::

      [ exclusions ]
      #include "exclusions.txt"

or copy and paste the content of that file to the exclusions section of
the gromacs topology file.

Statistical analysis
--------------------

For statistical analysis provides an interactive mode. To enter the
interactive mode, use the option followed by the file name of the
reference trajectory

::

      csg_boltzmann --top topol.tpr --trj traj.trr --cg mapping.xml

To get help on a specific command of the interactive mode, type

::

      help <command>

for example

::

      help hist
      help hist set periodic

Additionally, use the

::

      list

command for a list of available interactions. Note again that loads the
whole trajectory and all information on bonded interactions into the
memory. Hence, its main application should be single molecules. See the
introduction of this chapter for the command.

If a specific interaction shall be used, it can be referred to by

::

      molecule:interaction-group:index

Here, is the molecule number in the whole topology, is the name
specified in the section of the mapping file, and is the entry in the
list of interactions. For example, refers to the 10th bond named in
molecule 1. To specify a couple of interactions during analysis, either
give the interactions separated by a space or use wildcards (e.g. ).

To exit the interactive mode, use the command .

If analysis commands are to be read from a file, use the pipe or stdin
redirects from the shell.

::

      cat commands | csg_boltzmann topol.top --trj traj.trr --cg mapping.xml

Distribution functions and tabulated potentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Distribution functions (tabulated potentials) can be created with the
``hist`` (``tab``) command. For instance, to write out the distribution
function for all interactions of group AA-bond (where AA-bond is the
name specified in the mapping scheme) to the file AA.txt, type

::

      hist AA.txt *:AA-bond:*

The command

::

      hist set

prints a list of all parameters that can be changed for the histogram:
the number ``n`` of bins for the table, bounds ``min`` and ``max`` for
table values, scaling and normalizing, a flag ``periodic`` to ensure
periodic values in the table and an ``auto`` flag. If ``auto`` is set to
1, bounds are calculated automatically, otherwise they can be specified
by ``min`` and ``max``. Larger values in the table might extend those
bounds, specified by parameter ``extend``.

To directly write the Boltzmann-inverted potential, the ``tab`` command
can be used. Its usage and options are very similar to the ``hist``
command. If tabulated potentials are written, special care should be
taken to the parameters ``T`` (temperature) and the ``scale``. The
``scale`` enables volume normalization as given in . Possible values are
``no`` (no scaling), ``bond`` (normalize bonds) and ``angle`` (normalize
angles). To write out the tabulated potential for an angle potential at
a temperature of 300K, for instance, type:

::

      tab set T 300
      tab set scale angle
      tab angle.pot *:angle:*

The table is then written into the file in the format described in sec.
[sec:table\_formats]. An optional correlation analysis is described in
the next section. After the file has been created by command , the
potential is prepared for the coarse-grained run in chapter
[sec:usage:cgrun].

Correlation analysis
~~~~~~~~~~~~~~~~~~~~

The factorization of :math:`P` in assumed uncorrelated quantities.
offers two ways to evaluate correlations of interactions. One option is
to use the linear correlation coefficient (command ``cor``).

However, this is not a good measure since ``cor`` calculates the linear
correlation only which might often lead to misleading
results :raw-latex:`\cite{Ruehle:2009.a}`. An example for such a case
are the two correlated random variables :math:`X \sim U[-1,1]` with
uniform distribution, and :math:`Y:=X^2`. A simple calculation shows
:math:`cov(X,Y)=0` and therefore

.. math:: cor=\frac{cov(X,Y)}{\sqrt{var(X)var(Y)}}=0.

 A better way is to create 2D histograms. This can be done by specifying
all values (e.g. bond length, angle, dihedral value) using the command
*vals*, e.g.:

::

      vals vals.txt 1:AA-bond:1 1:AAA-angle:A

This will create a file which contains 3 columns, the first being the
time, and the second and third being bond and angle, respectively.
Columns 2 and 3 can either be used to generate the 2D histogram, or a
simpler plot of column 3 over 2, whose density of points reflect the
probability.

Two examples for 2D histograms are shown below: one for the propane
molecule and one for hexane.

.. figure:: fig/propane_hist2d
   :alt: hexane histograms: before and after the coarse-grained run
   :height: 2.80000cm

   hexane histograms: before and after the coarse-grained run

.. figure:: fig/hexane2
   :alt: hexane histograms: before and after the coarse-grained run
   :height: 3.60000cm

   hexane histograms: before and after the coarse-grained run

The two plots show the correlations between angle and bondlength for
both molecules. In the case of propane, the two quantities are not
correlated as shown by the centered distribution, while correlations
exist in the case of hexane. Moreover, it is visible from the hexane
plot that the partition of the correlations has changed slightly during
coarse-graining.

The tabulated potentials created in this section can be further modified
and prepared for the coarse-grained run: This includes fitting of a
smooth functional form, extrapolation and clipping of poorly sampled
regions. Further processing of the potential is decribed in chapter
[sec:usage:cgrun].

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

Iterative methods
=================

The following sections deal with the methods of Iterative Boltzmann
Inversion (), Inverse Monte Carlo (), and Relative Entropy ().

In general, , , and are implemented within the same framework.
Therefore, most settings and parameters of those methods are similar and
thus described in a general section (see sec.
[sec:iterative\_methods\_imc]). Further information on iterative methods
follows in the next chapters, in particular on the , , and methods.

.. figure:: usage/fig/flow_ibi.eps
   :alt: [fig:flow\_ibi]Flowchart to perform iterative Boltzmann
   inversion.
   :width: 7.00000cm

   [fig:flow\_ibi]Flowchart to perform iterative Boltzmann inversion.

Iterative workflow control
--------------------------

.. figure:: functionality/fig/flowchart.eps
   :alt: [fig:flowchart]Block-scheme of the workflow control for the
   iterative methods. The most time-consuming parts are marked in red.
   :width: 7.00000cm

   [fig:flowchart]Block-scheme of the workflow control for the iterative
   methods. The most time-consuming parts are marked in red.

Iterative workflow control is essential for the , , and methods.

The general idea of iterative workflow is sketched in
fig. [fig:flowchart]. During the global initialization the initial guess
for the coarse-grained potential is calculated from the reference
function or converted from a given potential guess into the internal
format. The actual iterative step starts with an iteration
initialization. It searches for possible checkpoints and copies and
converts files from the previous step and the base directory. Then, the
simulation run is prepared by converting potentials into the format
required by the external sampling program and the actual sampling is
performed.

After sampling the phasespace, the potential update is calculated.
Often, the update requires postprocessing, such as smoothing,
interpolation, extrapolation or fitting to an analytical form.

Finally, the new potential is determined and postprocessed. If the
iterative process continues, the next iterative step will start to
initialize.

How to start:
^^^^^^^^^^^^^

The first thing to do is generate reference distribution functions.
These might come from experiments or from atomistic simulations. To get
reasonable results out of the iterative process, the reference
distributions should be of good quality (little noise, etc).

can create initial guesses for the coarse-grained potentials by
boltzmann inverting the distribution function. If a custom initial guess
for an interaction shall be used instead, the table can be provided in
*:math:`<`\ interaction\ :math:`>`.pot.in*. As already mentioned,
automatically creates potential tables to run a simulation. However, it
does not know how to run a coarse-grained simulation. Therefore, all
files needed to run a coarse-grained simulation, except for the
potentials that are iteratively refined, must be provided and added to
the in the settings -file. If an atomistic topology and a mapping
definition are present, offers tools to assist the setup of a
coarse-grained topology (see chapter [sec:usage:cgrun]).

To get an overview of how input files look like, it is suggested to take
a look at one of the tutorials provided on .

In what follows we describe how to set up the iterative coarse-graining,
run the main script, continue the run, and add customized scripts.

Preparing the run
~~~~~~~~~~~~~~~~~

To start the first iteration, one has to prepare the input for the
sampling program. This means that all files for running a coarse-grained
simulation must be present and described in a separate file, in our case
``settings.xml`` (see sec. [sec:setting\_files] for details). An extract
from this file is given below. The only exception are tabulated
potentials, which will be created and updated by the script in the
course of the iterative process.

The input files include: target distributions, initial guess (optional)
and a list of interactions to be iteratively refined. As a target
distribution, any table file can be given (e.g. output from ``g_rdf``).
The program automatically takes care to resample the table to the
correct grid spacing according to the options provided in
``settings.xml``.

The initial guess is normally taken as a potential of mean force and is
generated by Boltzmann-inversion of the corresponding distribution
function. It is written in ``step_000/<name>.pot.new``. If you want to
manually specify the initial guess for a specific interaction, write the
potential table to a file called ``<name>.pot.in`` in the folder where
you plan to run the iterative procedure.

A list of interactions to be iteratively refined has to be given in the
options file. As an example, the ``setting.xml`` file for a propane is
shown in listing [list:settings]. For more details, see the full
description of all options in ref. [sec:ref\_options].

Starting the iterative process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After all input files have been set up, the run can be started by

::

      csg_inverse --options settings.xml

Each iteration is stored in a separate directory, named
``step_<iteration>``. ``step_000`` is a special folder which contains
the initial setup. For each new iteration, the files required to run the
CG simulation (as specified in the config file) are copied to the
current working directory. The updated potentials are copied from the
last step, ``step_<n-1>/<interaction>.pot.new``, and used as the new
working potentials ``step_<n>/<interaction>.pot.cur``.

After the run preparation, all potentials are converted into the format
of the sampling program and the simulation starts. Once the sampling has
finished, analysis programs generate new distributions, which are stored
in ``<interaction>.dist.new``, and new potential updates, stored in
``<interaction>.dpot.new``.

Before adding the update to the old potential, it can be processed in
the ``post_update`` step. For each script that is specified in the
postupdate, ``<interaction>.dpot.new`` is renamed to
``<interaction>.dpot.old`` and stored in
``<interaction>.dpot.<a-number>`` before the processing script is
called. Each processing script uses the current potential update
``<interaction>.dpot.cur`` and writes the processed update to
``<interaction>.dpot.new``. As an example, a pressure correction is
implemented as a postupdate script within this framework.

After all postupdate scripts have been called, the update is added to
the potential and the new potential ``<interaction>.pot.new`` is
written. Additional post-processing of the potential can be performed in
the ``post_add`` step which is analogous to the ``post_update`` step
except for a potential instead of an update.

To summarize, we list all standard output files for each iterative step:

+-----------------------+------------------------------------------------------------------------+
| ``*.dist.new``        | distribution functions of the current step                             |
+-----------------------+------------------------------------------------------------------------+
| ``*.dpot.new``        | the final potential update, created by ``calc_update``                 |
+-----------------------+------------------------------------------------------------------------+
| ``*.dpot.<number>``   | for each postupdate script, the ``.dpot.new`` is saved and a new one   |
+-----------------------+------------------------------------------------------------------------+
|                       | is created                                                             |
+-----------------------+------------------------------------------------------------------------+
| ``*.pot.cur``         | the current potential used for the actual run                          |
+-----------------------+------------------------------------------------------------------------+
| ``*.pot.new``         | the new potential after the add step                                   |
+-----------------------+------------------------------------------------------------------------+
| ``*.pot.<number>``    | same as ``dpot.<number>`` but for ``post_add``                         |
+-----------------------+------------------------------------------------------------------------+

If a sub-step fails during the iteration, additional information can be
found in the log file. The name of the log file is specified in the
steering file.

Restarting and continuing
~~~~~~~~~~~~~~~~~~~~~~~~~

The interrupted or finished iterative process can be restarted either by
extending a finished run or by restarting the interrupted run. When the
script is called, it automatically checks for a file called ``done`` in
the current directory. If this file is found, the program assumes that
the run is finished. To extend the run, simply increase in the settings
file and remove the file called ``done``. After that, can be restarted,
which will automatically recognize existing steps and continue after the
last one.

If the iteration was interrupted, the script might not be able to
restart on its own. In this case, the easiest solution is to delete the
last step and start again. The script will then repeat the last step and
continue. However, this method is not always practical since sampling
and analysis might be time-consuming and the run might have only crashed
due to some inadequate post processing option. To avoid repeating the
entire run, the script creates a file with restart points and labels
already completed steps such as simulation, analysis, etc. The file name
is specified in the option . If specific actions should be redone, one
can simply remove the corresponding lines from this file. Note that a
file ``done`` is also created in each folder for those steps which have
been successfully finished.

Iterative Boltzmann Inversion
-----------------------------

Input preparation
~~~~~~~~~~~~~~~~~

This section describes the usage of , implemented within the scripting
framework described in the previous section [sec:iterative\_workflow].
It is suggested to get a basic understanding of this framework before
proceeding.

An outline of the workflow for performing is given in .

To specify Iterative Boltzmann Inversion as algorithm in the script, add
``ibi`` in the ``method`` section of the setting file as shown below.

::

      <cg>
        ...
        <inverse>
          <method>ibi</method>
        </inverse>
      </cg>

Inverse Monte Carlo
-------------------

In this section, additional options are described to run coarse
graining. The usage of is similar to the one of and understanding the
use of the scripting framework described in
chapter [sec:iterative\_workflow] is necessary.

**WARNING: multicomponent is still experimental!**

General considerations
~~~~~~~~~~~~~~~~~~~~~~

In comparison to , needs significantly more statistics to calculate the
potential update:raw-latex:`\cite{Ruehle:2009.a}`. It is advisable to
perform smoothing on the potential update. Smoothing can be performed as
described in . In addition, can lead to problems related to finite size:
for methanol, an undersized system proved to lead to a linear shift in
the potential:raw-latex:`\cite{Ruehle:2009.a}`. It is therefore always
necessary to check that the system size is sufficiently large and that
runlength csg smoothing iterations are well balanced.

Correlation groups
~~~~~~~~~~~~~~~~~~

Unlike , also takes cross-correlations of interactions into account in
order to calculate the update. However, it might not always be
beneficial to evaluate cross-correlations of all pairs of interactions.
By specifying , allows to define groups of interactions, amongst which
cross-correlations are taken into account, where can be any name.

::

      <non-bonded>
        <name>CG-CG</name>
        <type1>CG</type1>
        <type2>CG</type2>
        ...
        <imc>
          <group>solvent</group>
       </imc>
      </non-bonded>
      <non-bonded>

Regularization
~~~~~~~~~~~~~~

To use the regularized version of IMC a :math:`\lambda` value :math:`>0`
has to be specified by setting . If set to :math:`0` (default value) the
unregularized version of IMC is applied.

::

     <non-bonded>
       <name>CG-CG</name>
       <type1>CG</type1>
       <type2>CG</type2>
        ...
       <inverse>
         <imc>
           <reg>300</reg>
         </imc>
       </inverse>
     </non-bonded>

Relative Entropy
----------------

In this section, additional options are described to run coarse
graining. The usage of is similar to the one of and and understanding
the use of the scripting framework described in
chapter [sec:iterative\_workflow] is necessary.

Currently, implementation supports optimization of two-body non-bonded
pair interactions. Support for bonded and N-body interactions is
possible by further extension of implementation.

Potential function and parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In , CG potentials are modeled using analytical functional forms.
Therefore, for each CG interaction, an analytical functional must be
specified in the setting file as

::

      <non-bonded>
        <name>CG-CG</name>
        <type1>CG</type1>
        <type2>CG</type2>
        ...
        <re>
          <function>cbspl or lj126</function>
            <cbspl>
              <nknots>48</nknots>
            </cbspl>
        </re>
        ...
      </non-bonded>

Currently, standard Lennard-Jones 12-6 (lj126) and uniform cubic
B-splines-based piecewise polynomial (cbspl) functional forms are
supported. For lj126, the parameters to optimize are the usual
:math:`C_{12}` and :math:`C_{6}`. The cbspl form is defined as

.. math::

   \label{eq:cbspl}
   u_{\text{cbspl}}(r) = \left[\begin{array}{cccc}
       1 & t & t^2 & t^3 \end{array}\right]
   \frac{1}{6}
   \left[ \begin{array}{rrrr}
       1 & 4 & 1 & 0 \\
       -3 & 0 & 3 & 0 \\
       3 & -6 & 3 & 0 \\
       -1 & 3 & -3 & 1 \end{array}\right]
   \left[ \begin{array}{l}
       c_{k} \\
       c_{k+1} \\
       c_{k+2} \\
       c_{k+3} \end{array}\right] ,

 where :math:`\{c_0,c_1,c_2,...,c_m\}` are the spline knot values
tabulated for :math:`m` evenly spaced intervals of size
:math:`\Delta r = r_{\text{cut}}/(m-2)` along the separation distance
:math:`r_{i} = i\times\Delta r` with the cut-off :math:`r_{\text{cut}}`,
and :math:`t` is given by

.. math::

   \label{eq:cbspl_t}
   t = \frac{r-r_{k}}{\Delta r} ,

 where index :math:`k` is determined such that
:math:`r_{k}\leq r < r_{k+1}`. For cbspl, the knot values,
:math:`\{c_0,c_1,c_2,...,c_m\}`, are optimized. The number of knot
values to use must be specified in the setting file as shown in the
above snippet. :math:`u_{\text{cbspl}}(r)` exhibits remarkable
flexibility, and it can represent various complex functional
characteristics of pair potentials for sufficiently large number of
knots.

Update scaling parameter
~~~~~~~~~~~~~~~~~~~~~~~~

Depending on the quality of the initial guess and sensitivity of the CG
system to the CG parameters, scaling of the parameter update size may be
required to ensure the stability and convergence of the minimization.
The scaling parameter, :math:`\chi\in(0...1)`, value can be specified in
the settings file.

Statistical averaging of parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Due to stochastic nature of the CG simulations, near convergence, the CG
potential paramters may fluctuate around the mean converged values.
Therefore, the optimal CG parameters can be estimated by averaging over
the last few iterations. To specify averaging, the ``average``, keyword
should be specified in the ``post_update`` options in the settings file.

General considerations
~~~~~~~~~~~~~~~~~~~~~~

To ensure the stability of the relative entropy minimization, some
precautionary measures are taken. For the Newton-Raphson update to
converge towards a minimum, the Hessian, :math:`\mathbf{H}`, must be
positive definite at each step. With a good initial guess for the CG
parameters and by adjusting the value of the relaxation parameter,
:math:`\chi`, stability of the Newton-Raphson method can be ensured. One
approach to initialize the CG parameters can be to fit them to PMF
obtained by inverting the pair distributions of the CG sites obtained
from the reference AA ensemble. For the lj126 and cbspl forms, which are
linear in its parameters, the second derivative of
:math:`S_{\text{rel}}` is never negative, hence the minimization
converges to a single global minimum. However, due to locality property
of the cbspl form, i.e., update to :math:`c_i` affects only the value of
the potential near :math:`r_i`, and the poor sampling of the very small
separation distances in the high repulsive core, the rows of
:math:`\mathbf{H}` corresponding to the first few spline knots in the
repulsive core may become zero causing :math:`\mathbf{H}` to be a
singular matrix. To avoid this singularity issue, we specify a minimum
separation distance, :math:`r_{\text{min}}`, for each CG pair
interaction and remove the spline knots corresponding to the
:math:`r\le r_{\text{min}}` region from the Newton-Raphson update. Once
the remaining knot values are updated, the knot values in the poorly
sampled region, i.e., :math:`r\le r_{\text{min}}`, are linearly
extrapolated. The value of :math:`r_{\text{min}}` at each iteration is
estimated from the minimum distance at which the CG RDF from the CG-MD
simulation is nonzero. Also, to ensure that the CG pair potentials and
forces go smoothly to zero near :math:`r_{\text{cut}}`, 2 knot values
before and after :math:`r_{\text{cut}}`, i.e., total 4, are fixed to
zero.

Pressure correction
-------------------

The pressure of the coarse-grained system usually does not match the
pressure of the full atomistic system. This is because iterative
Boltzmann inversion only targets structural properties but not
thermodynamic properties. In order correct the pressure in such a way
that it matches the target pressure ()., different strategies have been
used based on small modifications of the potential. The correction can
be enable by adding pressure to the list of scripts. The type of
pressure correction is selected by setting .

Simple pressure correction
~~~~~~~~~~~~~~~~~~~~~~~~~~

In ref.:raw-latex:`\cite{Reith:2003}` a simple linear attractive
potential was added to the coarse-grained potential

.. math:: \Delta V(r)=A \left( 1-\frac{r}{r_{cutoff}} \right) \,,

 with prefactor :math:`A`

.. math:: A = -\sign(\Delta P)0.1k_{B}T\min(1,|f\Delta P) \,,

 :math:`\Delta p=P_i-P_\text{target}`, and scaling factor :math:`f` and
:math:`P_\text{target}` can be specified in the settings file as and .

As an example for a block doing simple pressure correction, every third
interaction is

::

    <post_update>pressure</post_update>
    <post_update_options>
      <pressure>
        <type>simple</type>
        <do>0 0 1</do>
        <simple>
          <scale>0.0003</scale>
        </simple>
      </pressure
    </post_update_options>

Here, is the scaling factor :math:`f`. In order to get the correct
pressure it can become necessary to tune the scaling factor :math:`f`
during the iterative process.

Advanced pressure correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In :raw-latex:`\cite{Wang:2009}` a pressure correction based on the
virial expression of the pressure was introduced. The potential term
remains as in the simple form while a different sturcture of the
:math:`A` factor is used:

.. math:: A = \left[\frac{-2\pi\rho^{2}}{3r_{cut}}\int_{0}^{r_{cut}}r^{3}g_{i}(r)dr\right]A_{i}=\Delta P.

 This factor requires the particle density :math:` \rho ` as additional
input parameter, which is added as in the input file.

Kirkwood-Buff correction
------------------------

In order to reproduce the exact Kirkwood-Buff ingetrals (KBIs), an
correction term can be added into the coarse-grained
potential :raw-latex:`\cite{Ganguly:2012}`,

.. math:: \Delta U_{ij}^{(n)}(r) = \frac{k_{B}T}\;A\;(G_{ij}^{(n)} - G_{ij}^\text{ref})\left(1- \frac{r}{r_\text{ramp}}\right),

 where :math:`G_{ij}^{(ref)}` is the KBI calculated from the reference
all-atom simulation and :math:`G_{ij}^{(n)}` is the KBI after the
:math:`n^{th}` iteration.

The Kirkwood-Buff integrals are calculated from the radial distribution
functions as follows:

.. math::

   G_{ij} = 4\pi \int_0^\infty \left[ g_{ij}(r) - 1\right] r^2 dr~.
   \label{eq:kbi}

 For simulations of finite box size we calculate the running integral up
to distance :math:`R`

.. math:: G_{ij}(R) = 4\pi \int_0^R \left[ g_{ij}(r) - 1\right] r^2 dr~.

 The average of those running integrals in the interval, where
:math:`G_{ij}(R)` gets flat, gives a good estimate for :math:`G_{ij}`:

.. math:: G_{ij}\approx<G_{ij}(R)>|_{R=r_1}^{R=r_2}

 As an example for a block doing Kirkwood-Buff correction, every
iteraction without doing potential update

::

    <do_potential>0</do_potential>
    <post_update>kbibi</post_update>
    <post_update_options>
      <kbibi>
        <do>1</do>
        <start>1.0</start>
        <stop>1.4</stop>
        <factor>0.05</factor>
        <r_ramp>1.4</r_ramp>
      </kbibi>
    </post_update_options>

Here, is the scaling factor :math:`A`. is :math:`r_1` and is :math:`r_2`
used to calculate the average of :math:`G_{ij}(R)`.

Runtime optimization
--------------------

Most time per iteration is spent on running the coarse-grained system
and on calculating the statistics. To get a feeling on how much
statistics is needed, it is recommended to plot the distribution
functions and check whether they are sufficiently smooth. Bad statistics
lead to rough potential updates which might cause the iterative
refinement to fail. All runs should be long enough to produce
distributions/rdfs of reasonable quality.

Often, runtime can be improved by smoothing the potential updates. Our
experience has shown that it is better to smooth the potential update
instead of the rdf or potential itself. If the potential or rdf is
smoothed, sharp features like the first peak in water might get lost.
Smoothing on the delta potential works quite well, since the sharp
features are already present from the initial guess. By applying
iterations of a simple triangular smoothing
(:math:` \Delta U_i = 0.25 \Delta U_{i-1} + 0.5\Delta U_i + 0.25\Delta U_{i+1} `),
a reasonable coarse-grained potential for water could be produced in
less than 10 minutes. Smoothing is implemented as a post\_update script
and can be enabled by adding

::

      <post_update>smooth</post_update>
      <post_update_options>
        <smooth>
            <iterations>2</iterations>
        </smooth>
      </post_update_options>

to the inverse section of an interaction in the settings file.

Coordination Iterative Boltzmann Inversion
------------------------------------------

The method :math:`\mathcal{C}-`\ IBI (Coordination Iterative Boltzmann
Inversion) uses pair-wise cumulative coordination as a target function
within an iterative Boltzmann inversion. This method reproduces
solvation thermodynamics of binary and ternary mixtures
:raw-latex:`\cite{deOliveira:2016}`.

The estimation of coordination is given by:

.. math::

   \label{eq:coord}
   \mathcal{C}_{ij}(r) = 4\pi \int_{0}^{r} {\rm g}_{ij}(r')r'^{2}dr'

 with the indices :math:`i` and :math:`j` standing for every set of
pairs, uses a volume integral of :math:`{\rm g}(r)`.

The Kirkwood and Buff theory (KB) :raw-latex:`\cite{Kirkwood:1951}`
connects the pair-wise coordinations with particule fluctuations and,
thus, with the solution thermodynamics
:raw-latex:`\cite{Mukherji:2013,Naim:2006}`. This theory make use of the
Kirkwood-Buff integrals (KBI) :math:`{\rm G}_{ij}` defined as,

.. math::

   \label{eq:Gij}
   {\rm G}_{ij} = 4 \pi \int_{0}^{\infty} \left [ {\rm g}_{ij}(r) - 1 \right ] r^{2} dr.

 For big system sizes the :math:`{\rm G}_{ij}` can be approximated:

.. math::

   \label{eq:Gij_app}
   {\rm G}_{ij} = \mathcal{C}_{ij}(r) - \frac{4}{3} \pi r^{3},

 were the second therm is a volume correction to
:math:`\mathcal{C}_{ij}(r)`.

Thus the initial guess for the potential of the CG model is obtained
from the all atom simulations,

.. math::

   \label{eq:pot_ibi}
   {\rm V}_{0}(r) = -k_{B}T {\rm ln} \left [ {\rm g}_{ij}(r) \right ],

 however, the iterative protocol is modified to target
:math:`\mathcal{C}_{ij}(r)` given by,

.. math::

   \label{eq:pot_cibi}
   {\rm V}_{n}^{\mathcal{C}-{\rm IBI}}(r) = {\rm V}_{n-1}^{\mathcal{C}-{\rm IBI}}(r)
   + k_{B}T {\rm ln} \left [ \frac{\mathcal{C}_{ij}^{n-1}(r)}{\mathcal{C}_{ij}^{target}(r)} \right ].

To perform the :math:`\mathcal{C}-`\ IBI is necessary include some lines
inside of the .xml file:

::

     <cg>
      <non-bonded>
       <name>A-A</name>
       ...
       <inverse>
        <post_update>cibi</post_update>
        <post_update_options>
          <cibi>
            <do>1</do>
          </cibi>
        </post_update_options>
        ...
     </cg>

interface
=========

**WARNING: The interface is still experimental (in development) but it
does support the Iterative Boltzmann Inversion and Inverse Monte Carlo
schemes. The Force Matching might work as well, although it has not been
tested thoroughly.**

General remarks on using with 
------------------------------

The interface fully supports coarse-grain mapping of a full-atom system
previuosly simulated with any version of , including -Classic. However,
the full optimization of the effective potentials with the aid of
iterative methods will only become possible when the new release of
(4.06) is made public; the reason being the incapability of earlier
versions of using user-specified tabulated force-fields for
intramolecular, aka “bonded”, interactions: bonds, angles, dihedral
angles (torsions). Below the coarse-graining and CG force-field
optimization with the aid of the latest version (4.06+) are outlined.

Running with as MD simulation engine is very similar to doing so with .
The three types of required input files in the case of are: CONTROL –
containing the simulation directives and parameters (instead of ``.mdp``
file for ), FIELD – the topology and force-field specifications (instead
of ``.top`` and ``.tpr`` files), and CONFIG (instead of ``.gro`` file) –
the initial configuration file, containing the MD cell matrix and
particle coordinates (it can also include initial velocities and/or
forces); for details see manual. Most of the tools and scripts described
above in the case of using will work in the same manner, with the
following conventional substitutions for the (default) file names used
in options for scripts, as necessary:

::

    .dlpf = the topology read from FIELD or written to FIELD_CGV
    .dlpc = the configuration read from CONFIG or written to CONFIG_CGV
    .dlph = the trajectory read from HISTORY or written to HISTORY_CGV

It is also possible to specify file names different from the standard
convention, in which case the user has to use the corresponding
dot-preceded extension(s); for example: FA-FIELD.dlpf instead of FIELD
or CG-HISTORY.dlph instead of HISTORY\_CGV (see section
[sec:ref\_programs], as well as the man pages or output of commands,
with option ).

follows the conventions for file names and formats. Thus,
``csg_dlptopol`` and ``csg_map`` produce the CG topology (FIELD\_CGV by
default), configuration (CONFIG\_CGV), and/or trajectory (HISTORY\_CGV)
files fully compatible with and usable by . **Note that the ability of
these tools to read and write a plethora of different file formats
provides means to convert input and output files between the simulation
packages supported by , e.g. – or vice versa. The user is, however,
strongly advised to check the resulting files for consistency before
using them).**

Similarly, the distribution analysis and potential/force generation
utilities, such as ``csg_stat`` and scripts, will read and write
-formatted files; in particular, the tabulated force-field files
containing the potential and force/virial data: TABLE – for short-range
(VdW) “non-bonded” interactions, TABBND, TABANG and TABDIH – for
“bonded” interations: bonds, bending angles and dihedrals,
correspondingly (for the format details see manual). Note, however, that
the latter three files can only be used by (4.06+).

The user is advised to search for “dlpoly” through the
``csg_defaults.xml``, ``csg_table`` files and in scripts located in
``share/votca/scripts/inverse/`` in order to find out about the xml-tags
and options specific for ; see also sections [sec:ref\_options]
and [sec:csg\_table].

.. |image| image:: usage/fig/flow_boltzmann.eps
   :width: 7.00000cm
