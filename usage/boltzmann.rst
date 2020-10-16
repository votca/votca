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

.. |image| image:: usage/fig/flow_boltzmann.eps
   :width: 7.00000cm
