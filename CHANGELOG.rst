For more detailed information about the changes see the history of the
`repository <https://github.com/votca/csg/commits/master>`__.

Version 2021-rc.2 (released XX.01.21)
=====================================

Version 2021-rc.1 (released 15.01.21)
=====================================

-  remove topology item (#486)
-  switch to local includes (#487)
-  add default units to readers (#483)
-  update bead related classes (#475, #474, #469, #468, #460)
-  update boundary class (#457)
-  fix usage of PROJECT\_SOURCE\_DIR in cmake (#517)
-  document cmake options and packages (#520)
-  drop WITH\_GMX in cmake in favor of auto-detect (#523)
-  standardising include order and style in header files (#531)
-  make scripts always executable in builddir (#533, #550, #556)
-  rename REGRESSIONTEST\_TOLERANCE to INTEGRATIONTEST\_TOLERANCE (#534)
-  standardising include order and style in library files (#542)
-  move CI to GitHub Ations (#551, #552, #553, #561, #562, #593, #617,
   #632)
-  format CMake internal cpp code (#557)
-  update codacy badge (#558)
-  convert markdown to rst format (#559)
-  read units stored in H5MD file (#549)
-  enable CXX only in CMake (#560)
-  refactored tests to use files instead of data in source code (#572)
-  corrected max index in csg_resample (#575)
-  add periodic extrapolation for dihedrals (#576)
-  bump required CMake version to 3.12 (#599)
-  Fixed boost test deprecation warnings (#598)
-  improve error message for broken binaries (#601)
-  Fixed memory leaks in tests and neighbourlist (#602)
-  format with clang-11 (#604)
-  use unique_ptr in factory (#603)
-  Convert remaining md files to rst (#616)
-  merge csg-manual (#621, #622)
-  merge in csgapps repo (#623)
-  Pdbs with more than 10000 resids can now be written (#627)

Version 1.6.4 (released 12.01.21)
=================================
* fix exclusion creation in csg_boltzmann (#624, #625)

Version 1.6.3 (released 09.12.20)
=================================

-  fix test dependencies for parallel ctest (#585)
-  fix trigger for gmx-2020 warning in run_gromacs.sh (#591)
-  switch to ghcr.io for CI (#595)
-  format code with clang-11 (#605)
-  add gmx2021 builds to CI (#607)
-  remove set-env call from CI (#608)
-  use master .clang-format in format action (#610)
-  strip windows line endings in readers (#614)

Version 1.6.2 *SuperGitta* (released 22.08.20)
=================================

-  move CI to GitHub Actions (#563, #567, #569)
-  add max scaling parameter to pressure correction (#565)
-  fix test dependencies for parallel ctest (#573)
-  fix multicomponent imc (#579, #581)

Version 1.6.1 (released 21.06.20)
=================================

-  add support for cma-3 (#158)
-  allow special characters in build directory (#521)
-  fix CI on Ubuntu-20.04 (#526)
-  drop support for gmx 2016 (#529, #547)
-  fix csg\_reupdate and add tests (#525)
-  fix detection of lmp and gmx in cmake (#540)
-  add test for lammpsdata reader (#544)

Version 1.6 *SuperPelagia* (released 17.04.20)
==============================================

-  enforce periodic boundaries for dihedrals (#500)
-  add warning about dihedrals (#500)
-  fix typo in option in potential\_extrapolate (#502)
-  update gmx tab. potential warning (#504)
-  fix exit code with clean for potential\_to\_lammps.sh (#506)
-  remove obsolete link to MKL from cmake (#508)
-  fix build with Gromacs 2021 (#510, #511)
-  warn about too long interaction in csg\_stat (#512)
-  fix testing in postadd\_plot.sh (#515)

Version 1.6\_rc2 (released 10.02.20)
====================================

-  fix inheritance bug in cmake files (#466)
-  fix build on 32-bit archs (#473)
-  fixed alignment problem in potential class (#471)
-  fix h5md reader (#478, #481)
-  added more tests for gmx (#489) and gmx-2020 support (#488)
-  add support for txt2tags-3.4 (#491)
-  fix cma seed for testing (#494)
-  remove old doxygen target (#498)

Version 1.6\_rc1 (released 04.12.19)
====================================

-  Major CMake refactor
-  Improved file readers
-  Added more tests
-  Added csg\_imc\_solve, removing octave and python solvers

Version 1.5.1 (released 20.11.19)
=================================

-  bug fixes for gcc-5
-  some bugfixes for the lammps reader
-  support dynamics box size in h5md
-  fix build on CentOs7

Version 1.5 *SuperVictor* (released 31.01.19)
=============================================

-  enable gitlab CI

Version 1.5\_rc3 (released 19.01.19)
====================================

-  travis: fixed bug in building tags
-  cmake: only run csg\_dump test if bash was found

Version 1.5\_rc2 (released 16.01.19)
====================================

-  added REGRESSIONTEST\_TOLERANCE cmake option
-  added bug fix for memory leak in bead class
-  added topology unit test

Version 1.5\_rc1 (released 28.12.18)
====================================

-  added support for 3 body distribution
-  allow usage of \_ in gmx options
-  support espresso-4.0
-  added LAMMPS data topology reader
-  added LAMMPS dump writer
-  added support for gmx-2019
-  move from gsl to Eigen3
-  added support for running gmx in multidir mode
-  extended pdb reader to read CONECT section
-  reworked perl interface to fix #179
-  make most install locations configureable
-  allow to change pot extrapolation fct for LAMMPS
-  dropped support for gmx-5.1
-  added unit and regression tests

Version 1.4.1 (released 02.09.17)
=================================

-  fix pkg-config files
-  fix build with gmx-2017
-  added CSG\_MDRUN\_OPTS to run\_gromacs.sh
-  table\_to\_tab.pl: fix forces for LAMMPS
-  csg\_inverse: warn users of bug #179
-  run\_gromacs.sh: check nstxout-compressed
-  tools: fix parallel build with mkl

Version 1.4 *SuperKurt* (released 29.10.16)
===========================================

-  added cibi tutorial and manual section
-  install mkl headers if mkl is enabled
-  updated copyright
-  xmltopologyparser: tokenize on \\n and \\t as well (#195)
-  tools: added support for boost-1.62

Version 1.4\_rc1 (released 26.09.16)
====================================

-  switch to c++-11
-  dihedral support in csg\_fmatch
-  support for tabulated bond in >=gmx-5
-  added full featured XMLTopologyReader
-  added regularization for IMC
-  changed neighborlist from std::list to std::vector
-  added cibi method
-  tools: added support for mkl as gsl replacement
-  lots of reader bug fixes
-  dropped support for gromacs-4 and clean up
-  dropped multi\_g\_rdf script
-  dropped thermforce iteration method
-  moved h5md reader to hdf5 without c++ api

Version 1.3.1 (released 19.08.16)
=================================

-  histogram: use floor() on bin value
-  calculator: fixed namespace
-  VOTCARC: added shebang
-  fixed gromacs detection with >=cmake-3.4

Version 1.3 *SuperUzma* (released 15.01.16)
===========================================

-  re-implemented csg\_boltzmann --excl
-  added support for upcoming gromacs 2016

Version 1.3\_rc1 (released 23.09.15)
====================================

-  added new iterative methods: relative entropy, simplex optimization
-  added support for using the following with iterative methods:
   hoomd-blue, lammps, ESPResSo, ESPResSo++, dl\_poly
-  added pre-simulation feature for GROMACS (e.g. for minimization)
-  added rudimentary support for IBI with bonded interaction
-  made pdb reader work with libgmx
-  added support for h5md, dl\_ploy file format
-  added support for numpy in IMC
-  cmake: added BUILD\_MANPAGES option, git support minor fixes
-  cmake: dropped internal boost replacement
-  many many many small bug fixes and improvements

Version 1.2.4 (released 31.08.14)
=================================

-  support for Gromacs 5.0
-  support for Boost 1.53
-  fixed use of nawk instead of gawk under MacOs
-  fixed python shebang
-  fixed linking issue under Fedora
-  fixed thermforce calculation for xsplit case

Version 1.2.3 (released 14.08.12)
=================================

-  improved AIX support
-  fixed install on 64-bit linux systems
-  fixed a bug in histogram class
-  fixed rdf calculation for r\_min > 0 (histogram bug)
-  updated documentation

Version 1.2.2 (released 10.01.12)
=================================

-  added numpy solver for IMC
-  cmake: updated FindGROMACS.cmake
-  fixed coredump in csg\_property (issue 114)
-  fixed namespace in Fedora
-  fixed problem with newlines in csg\_property
-  cmake: allow static fftw and gsl
-  added dummy c function for cmake
-  fixed conflicting type headers (real was defined)

Version 1.2.1 (released 25.08.11)
=================================

-  csg\_inverse: improve initial guess of the potential
-  csg\_inverse: fixes for min!=0
-  table\_extrapolate.pl: fixed flags and first point
-  fixed tf iteration for multiple components
-  fixed round-off error in grid search and csg\_calc
-  csg\_inverse: typo fixed and additional checks
-  fixed soname of libs
-  improved cmake checks and error messages
-  fixed pkg-config file

Version 1.2 *SuperDoris* (released 17.06.11)
============================================

-  changed buildsystem to cmake
-  added thermforce iteration method
-  added csg\_density
-  a lot of framework clean up
-  added type selector name:\*
-  allow long and restart of simulations
-  added database class through sqlite3

Version 1.1.2 (released 04.04.11)
=================================

-  csg\_fmatch: added support for known forces (--trj-force option)
-  fixed head of votca.7 manpage

Version 1.1.1 (released 01.03.11)
=================================

-  fixed csg\_inverse --clean
-  make postupdate pressure work again
-  fixed bug when reading exclusions from tpr
-  end with error in csg\_stat if bead type does not exist (issue 77)

Version 1.1 *SuperAnn* (released 18.02.11)
==========================================

-  added support for gromacs 5.0
-  csg\_dump: can dump exclusion
-  added boundarycondition class
-  added man pages, man7 and man1 for all bins
-  csg\_inverse: renamed ibm to ibi
-  csg\_inverse: many internal improvements
-  csg\_stat: added: thread support, read exclusions from tpr file, uses
   grid search by default
-  csg\_inverse: added: convergence check, postadd plot, better logging,
   weaker die
-  csg\_resample: added boundary conditions option and akima spline
   support
-  csg\_stat or csg\_fmatch give an error if trj not given (issue 29)
-  csg\_get\_interaction\_property knows about defaults
-  fixed segfault in mapping (Fixes issue 27)
-  fixed bug in gromacs writer (frame.bX = true)
-  fixed segfault in gromacs writer (issue 54)
-  added thread class
-  added spline class, with akima spline, linear spline
-  random.cc: avoid calling of exit()
-  added lexical cast class

Version 1.0.1 (released 01.12.10)
=================================

-  fixed custom md programs in sim scripts (issue 1)
-  completion file is back from tools
-  issue #21: fixed strange kink when pot.in was provided
-  added --disable-rc-files to configure
-  csg\_call/csg\_inverse: added installdir as failback for CSGSHARE
-  fixed a bug in VOTCARC.csh for empty LD\_LIBRARY\_PATH
-  completion file has moved back to csg
-  added --disable-rc-files to configure
-  updated bundled libtool to 2.2.10

Version 1.0 (released 30.09.10)
===============================

-  added postupdate script for scaling the update
-  imc and csg\_stat: no longer require dummy mapping file (--no-map)
   option
-  allow comments in tables
-  fixed bug in pressure correction when p is negative
-  added support for gromacs devel version
-  fixed a bug when compiling with gcc-4.4
-  fixed a bug that pot.cur was change at every step
-  added application class for easy implementation of analysis programs
-  fixed bug if initial potential was given and not used
-  restart points are no longer deleted after step finished
-  csg\_inverse: preliminary reader for ESPResSo Blockfiles and ESPResSo
-  preliminary reader for LAMMPS dump files (very limited features)
-  allow compling without gromacs
-  a lot new xml optionsfull support for gromacs 4.5
-  added libvotca\_expat to allow compiling without expat
-  allow comments in tables
-  added application class to create standardized applications
-  all boost dependecy are now in tools
-  fixes in table format, flags is always last row now
-  allow compling without fftw (needed for csg\_boltzmann only)
-  allow compling without gsl (needed for csg\_resample and csg\_fmatch)

Version 1.0\_rc5 (released 16.03.10)
====================================

-  fixed --first-frame option (--first-frame 1 before started at second
   frame)
-  fixed compatibility issue when using gromacs development version
-  updated configure, see --help
-  added multi\_g\_density
-  CSGRC is replaced by VOTCARC of votca\_tools
-  using libexpat instead of libxml2
-  added libvotca\_boost to allow compiling without boost
-  using pkg-config to detect package flags
-  compiles under AIX with xlC
-  added VOTCARC to initialize all votca parts
-  updated configure, see --help

Version 1.0\_rc4 (released 08.02.10)
====================================

-  using libtool to build shared libs\\
-  fixed a bug in error calculation of multi\_g\_rdf

Version 1.0\_rc3 (released 29.01.10)
====================================

-  added option --wall-time to csg\_inverse if run on a queueing system
-  added option for IBI to run in parallel
-  multi\_g\_rdf, a multiplexed version of g\_rdf was added
-  added some options to csg\_call
-  csg\_resample now also calc derivatives
-  fixed a bug in reading stuff from mpd file
-  corrected bug in tokenizer
-  fixed a bug in calculation of version string
-  some fixes concerning autotools

Version 1.0\_rc2 (released 16.12.09)
====================================

-  added version string to scripts
-  fixed typo in calculation of version string
-  added NOTICE and LICENSE to the dist tarball

Version 1.0\_rc1 (released 11.12.09)
====================================

-  initial version
