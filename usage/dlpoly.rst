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
