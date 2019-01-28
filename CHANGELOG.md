For more detailed information about the changes see the history of the [repository](https://github.com/votca/csg/commits/master).

## Version 1.5 (released XX.01.19)
 * enable gitlab CI

## Version 1.5_rc3 (released 19.01.19)
 * travis: fixed bug in building tags
 * cmake: only run csg_dump test if bash was found

## Version 1.5_rc2 (released 16.01.19)
* added REGRESSIONTEST_TOLERANCE cmake option
* added bug fix for memory leak in bead class
* added topology unit test

## Version 1.5_rc1 (released 28.12.18)
* added support for 3 body distribution
* allow usage of _ in gmx options
* support espresso-4.0
* added LAMMPS data topology reader
* added LAMMPS dump writer
* added support for gmx-2019
* move from gsl to Eigen3
* added support for running gmx in multidir mode
* extended pdb reader to read CONECT section
* reworked perl interface to fix #179
* make most install locations configureable
* allow to change pot extrapolation fct for LAMMPS
* dropped support for gmx-5.1
* added unit and regression tests

## Version 1.4.1 (released 02.09.17)
* fix pkg-config files
* fix build with gmx-2017
* added CSG_MDRUN_OPTS to run_gromacs.sh
* table_to_tab.pl: fix forces for LAMMPS
* csg_inverse: warn users of bug #179
* run_gromacs.sh: check nstxout-compressed
* tools: fix parallel build with mkl

## Version 1.4 _SuperKurt_ (released 29.10.16)
* added cibi tutorial and manual section
* install mkl headers if mkl is enabled
* updated copyright
* xmltopologyparser: tokenize on \n and \t as well (#195)
* tools: added support for boost-1.62

## Version 1.4_rc1 (released 26.09.16)
 * switch to c++-11
 * dihedral support in csg_fmatch
 * support for tabulated bond in >=gmx-5
 * added full featured XMLTopologyReader
 * added regularization for IMC
 * changed neighborlist from std::list to std::vector
 * added cibi method
 * tools: added support for mkl as gsl replacement
 * lots of reader bug fixes
 * dropped support for gromacs-4 and clean up
 * dropped multi_g_rdf script
 * dropped thermforce iteration method
 * moved h5md reader to hdf5 without c++ api

## Version 1.3.1 (released 19.08.16)
* histogram: use floor() on bin value
* calculator: fixed namespace
* VOTCARC: added shebang
* fixed gromacs detection with >=cmake-3.4

## Version 1.3 _SuperUzma_ (released 15.01.16)
* re-implemented csg_boltzmann --excl
* added support for upcoming gromacs 2016

## Version 1.3_rc1 (released 23.09.15)

* added new iterative methods: relative entropy, simplex optimization
* added support for using hoomd-blue, lammps, ESPResSo, ESPResSo++ and dl_poly for iterative methods
* added pre-simulation feature for GROMACS (e.g. for minimization)
* added rudimentary support for IBI with bonded interaction
* made pdb reader work with libgmx
* added support for h5md, dl_ploy file format
* added support for numpy in IMC
* cmake: added BUILD_MANPAGES option, git support minor fixes
* cmake: dropped internal boost replacement
* many many many small bug fixes and improvements

## Version 1.2.4 (released 31.08.14)

* support for Gromacs 5.0
* support for Boost 1.53
* fixed use of nawk instead of gawk under MacOs
* fixed python shebang
* fixed linking issue under Fedora
* fixed thermforce calculation for xsplit case

## Version 1.2.3 (released 14.08.12)

* improved AIX support
* fixed install on 64-bit linux systems
* fixed a bug in histogram class
* fixed rdf calculation for r_min > 0 (histogram bug)
* updated documentation

## Version 1.2.2 (released 10.01.12)

* added numpy solver for IMC
* cmake: updated FindGROMACS.cmake
* fixed coredump in csg_property (issue 114)
* fixed namespace in Fedora
* fixed problem with newlines in csg_property
* cmake: allow static fftw and gsl
* added dummy c function for cmake
* fixed conflicting type headers (real was defined)

## Version 1.2.1 (released 25.08.11)

* csg_inverse: improve initial guess of the potential
* csg_inverse: fixes for min!=0
* table_extrapolate.pl: fixed flags and first point
* fixed tf iteration for multiple components
* fixed round-off error in grid search and csg_calc
* csg_inverse: typo fixed and additional checks
* fixed soname of libs
* improved cmake checks and error messages
* fixed pkg-config file

## Version 1.2 _SuperDoris_ (released 17.06.11)

* changed buildsystem to cmake
* added thermforce iteration method
* added csg_density
* a lot of framework clean up
* added type selector name:*
* allow long and restart of simulations
* added database class through sqlite3

## Version 1.1.2 (released 04.04.11)

* csg_fmatch: added support for known forces (--trj-force option)
* fixed head of votca.7 manpage

## Version 1.1.1 (released 01.03.11)

* fixed csg_inverse --clean
* make postupdate pressure work again
* fixed bug when reading exclusions from tpr
* end with error in csg_stat if bead type does not exist (issue 77)

## Version 1.1 _SuperAnn_ (released 18.02.11)

* added support for gromacs 5.0
* csg_dump: can dump exclusion
* added boundarycondition class
* added man pages, man7 and man1 for all bins
* csg_inverse: renamed ibm to ibi
* csg_inverse: many internal improvements
* csg_stat: added thread support, read exclusions from tpr file, uses grid search by default
* csg_inverse: added convergence check, postadd plot, better logging and weaker die
* csg_resample: added boundary conditions option and akima spline support
* csg_stat or csg_fmatch give an error if trj not given (issue 29)
* csg_get_interaction_property knows about defaults
* fixed segfault in mapping (Fixes issue 27)
* fixed bug in gromacs writer (frame.bX = true)
* fixed segfault in gromacs writer (issue 54)
* added thread class
* added spline class, with akima spline, linear spline
* random.cc: avoid calling of exit()
* added lexical cast class

## Version 1.0.1 (released 01.12.10)

* fixed custom md programs in sim scripts (issue 1)
* completion file is back from tools
* issue #21: fixed strange kink when pot.in was provided
* added --disable-rc-files to configure
* csg_call/csg_inverse: added installdir as failback for CSGSHARE
* fixed a bug in VOTCARC.csh for empty LD_LIBRARY_PATH
* completion file has moved back to csg
* added --disable-rc-files to configure
* updated bundled libtool to 2.2.10

## Version 1.0 (released 30.09.10)

* added postupdate script for scaling the update
* added options to perform analysis without giving a mapping file (--no-map). No need for dummy mapping file in IMC and csg_stat anymore.
* allow comments in tables
* fixed bug in pressure correction when p is negative
* added support for gromacs devel version
* fixed a bug when compiling with gcc-4.4
* fixed a bug that pot.cur was change at every step
* added application class for easy implementation of analysis programs
* fixed bug if initial potential was given and not used
* restart points in iterative procedure are not anymore deleted after step finished
* preliminary reader for ESPResSo Blockfiles and ESPResSo as Simulation program for csg_inverse
* preliminary reader for LAMMPS dump files (very limited features)
* allow compling without gromacs
* a lot new xml optionsfull support for gromacs 4.5
* added libvotca_expat to allow compiling without expat
* allow comments in tables
* added application class to create standardized applications
* all boost dependecy are now in tools
* fixes in table format, flags is always last row now
* allow compling without fftw (needed for csg_boltzmann only)
* allow compling without gsl (needed for csg_resample and csg_fmatch)

## Version 1.0_rc5 (released 16.03.10)

* fixed --first-frame option (--first-frame 1 before started at second frame)
* fixed compatibility issue when using gromacs development version
* updated configure, see --help
* added multi_g_density
* CSGRC is replaced by VOTCARC of votca_tools
* using libexpat instead of libxml2
* added libvotca_boost to allow compiling without boost
* using pkg-config to detect package flags
* compiles under AIX with xlC
* added VOTCARC to initialize all votca parts
* updated configure, see --help

## Version 1.0_rc4 (released 08.02.10)

* using libtool to build shared libs\
* fixed a bug in error calculation of multi_g_rdf

## Version 1.0_rc3 (released 29.01.10)

* added option --wall-time to csg_inverse if run on a queueing system
* added option for IBI to run in parallel
* multi_g_rdf, a multiplexed version of g_rdf was added
* added some options to csg_call
* csg_resample now also calc derivatives
* fixed a bug in reading stuff from mpd file
* corrected bug in tokenizer
* fixed a bug in calculation of version string
* some fixes concerning autotools

## Version 1.0_rc2 (released 16.12.09)

* added version string to scripts
* fixed typo in calculation of version string
* added NOTICE and LICENSE to the dist tarball

## Version 1.0_rc1 (released 11.12.09)

* initial version

