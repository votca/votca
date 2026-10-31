Advanced topics
===============

Used external packages
----------------------

GROMACS
~~~~~~~

`Link to Homepage <https://manual.gromacs.org/documentation/>`_

DL\_POLY
~~~~~~~~

`Link to Homepage <https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx>`_


DL\_POLY interface
------------------

**WARNING: The DL\_POLY interface is experimental and has not been tested
thoroughly.**

The DL\_POLY interface supports coarse-grain mapping of a full-atom
system previously simulated with any version of DL\_POLY, including
DL\_POLY-Classic. The three types of DL\_POLY input files are: CONTROL –
containing the simulation directives and parameters (instead of ``.mdp``
file for GROMACS), FIELD – the topology and force-field specifications
(instead of ``.top`` and ``.tpr`` files), and CONFIG (instead of ``.gro``
file) – the initial configuration file, containing the MD cell matrix and
particle coordinates (it can also include initial velocities and/or
forces); for details see the DL\_POLY-4 manual. VOTCA tools recognize
DL\_POLY files by the following extensions:

.. code:: none

    .dlpf = the topology read from FIELD or written to FIELD_CGV
    .dlpc = the configuration read from CONFIG or written to CONFIG_CGV
    .dlph = the trajectory read from HISTORY or written to HISTORY_CGV

It is also possible to specify file names different from the standard
DL\_POLY convention, in which case the user has to use the corresponding
dot-preceded extension(s); for example: FA-FIELD.dlpf instead of FIELD
or CG-HISTORY.dlph instead of HISTORY\_CGV (see
:ref:`reference_programs`, as well as the man pages or output of
VOTCA commands, with option ``--help``).

VOTCA follows the DL\_POLY conventions for file names and formats. Thus,
``csg_dlptopol`` and ``csg_map`` produce the CG topology (FIELD\_CGV by
default), configuration (CONFIG\_CGV), and/or trajectory (HISTORY\_CGV)
files fully compatible with and usable by DL\_POLY. **Note that the
ability of these tools to read and write a plethora of different file
formats provides means to convert input and output files between the
simulation packages supported by VOTCA, e.g. GROMACS – DL\_POLY or vice
versa. The user is, however, strongly advised to check the resulting
files for consistency before using them.**

VOTCA does not write DL\_POLY tabulated force-field files (TABLE, TABBND,
TABANG, TABDIH); tables from ``csg_stat`` or ``csg_fmatch`` have to be
converted with other tools.
