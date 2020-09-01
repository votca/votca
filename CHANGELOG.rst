For more detailed information about the changes see the history of the
`repository <https://github.com/votca/xtp/commits/master>`__.

Version 2021-dev
===============

-  G0W0 rpa energies are fixed in BSE (#383)
-  Hqp part of BSE can be separate from QP range (#388)
-  Fixed Cuda implementation (#391)
-  Plain, linear, and Anderson mixing in evGW
-  Use a general syntax to represent all QMPackages input (#318)
-  Remove support for both Gaussian and NWChem (#318)
-  Fixed executable path check (#400)
-  Usage of offdiagonal elements of Hqp in BSE optional, default: with
   offdiagonals (#402)
-  refactored MO reordering of external QMPackages (#406)
-  Add defaults to XML (#319, #413)
-  fix build with gcc-10 (#422, #424)
-  Add name command line argument (#426)
-  fix usage of PROJECT\_SOURCE\_DIR in cmake (#431)
-  check that temperature is positive for KMC (#432)
-  document cmake options and packages (#434)
-  remove long gone MKL option for gitlab-ci (#435)
-  add documentation using Sphinx (#433)
-  add kokkos-based Ewald solver (#444, #446, #485)
-  allow multiple choices for the calculator options (#442, #445)
-  move bsecoupling defaults to xml (#451)
-  address missing includes (#452)
-  fix readjobfiles path (#453)
-  standardising include order and style in header files (#448)
-  changed to OpenMP reductions instead of hand crafted solutions (#466,
   #471, #493)
-  switch to GitHub Actions as CI (#467, #474, #477, #482, #500, #502,
   #518, #549)
-  move tutorial to its own repo (#450)
-  removed share/packages folder (#475)
-  copy orca \*.gbw files in eqm (#472)
-  improve mkl support (#479)
-  add def2-tzvp, def2-tzvpp, def2-qzvp, def2-qzvpp for Rb,I,Ag,Xe
   (#480)
-  create a map between orca and libxc functional names (#461)
-  fix path to share/data (#486, #487)
-  remove special XML symbols (#488)
-  standardize includes in src files (#462)
-  add cmake exported target (#489)
-  update codacy badge (#490)
-  add mkl builds to CI (#496)
-  convert markdown to rst (#497)
-  enable CXX only in CMake (#499)
-  implement RPA correlation energy (#498)
-  added verbose option for rootfinder (#503)
-  fix orca pointcharges (#504)
-  remove constant cutoff from neighborlist (#508, #509)
-  fix default and flags properties issue (#511)
-  add more checks to the dft_parse_part (#510)
-  add NOSCF option for iqm in internal and external DFT (#483)
-  added versions to output files (#523)
-  write an empty state file if it does not exist yet (#526)
-  fix double unit conversion (#531)
-  Unittests can read input and reference values from files (#529)
-  fix double SCF in ORCA keyword line (#533)
-  xtp internally uses only single shells. (#539)
-  Refactor orbital reordering (#537)
-  fix build with clang-11 (#546)
-  Add a molden parser to xtp (#547)
   
Version 1.6.2 *SuperGitta* (released 22.08.20)
=================================

-  move CI to GitHub Actions (#512, #514, #516,
   #519)
-  add more checks to the dft_parse_part (#521)
-  make unit tests work in parallel (#530)

Version 1.6.1 (released 21.06.20)
=================================

-  fix warnings on Ubuntu 20.04 (#438, #460)
-  remove obsolete boost define (#440)

Version 1.6 *SuperPelagia* (released 17.04.20)
==============================================

-  fix 32-bit build (#381, #380)
-  remove duplicated basissets (#384, #386, #387)
-  fix clang-10 warnings (#394)
-  fix unit\_test\_eeinteractor on OpenSUSE (#341, #428)

Version 1.6\_rc2 (released 10.02.20)
====================================

-  fix remove giant logo from tarball (#337)
-  fix assertions related to GLIBCXX\_ASSERTIONS (#345)
-  remove unused boost serialisation (#346)
-  fix build on 32-bit archs (#347)
-  add ENABLE\_HIGH\_MEMORY\_TESTS cmake option (#356)
-  fix copyright (#363)
-  remove old doxygen target (#365)
-  fix some gcc10 warnings (#376)
-  Add external fields to dft calculations (#351, #353)
-  added def2 and cc basis sets (#355)
-  added apdft (#350)
-  added test to cubefile reader (#344)
-  fix state tracker (#333)
-  grid class refator (#335)
-  changed ppm screening (#371)

Version 1.6\_rc1 (released 04.12.19)
====================================

-  completely new statefile in hdf5 format for larger systems
-  new electrostatics with PCG solver
-  new QM/MM engine with freely configurable regions
-  exact GW-BSE for small systems
-  new iterative matrix solvers for large systems
-  CUDA support for parts of GW-BSE
-  full LAMMPS support
-  improved testing
-  new tutorial with more functionality
-  deleted netbeans support
-  removed old manual
-  Major CMake refactor

Version 1.5.1 (released 20.11.19)
=================================

-  remove exit() calls in the library
-  fix build on CentOs7

Version 1.5 *SuperVictor* (released 31.01.19)
=============================================

-  enable gitlab CI

Version 1.5\_rc3 (released 19.01.19)
====================================

-  travis: fixed bug in building tags

Version 1.5\_rc2 (released 16.01.19)
====================================

-  fix parallel build of manual
-  fix usage on inkscape on arm arch
-  clean up namespace in header

Version 1.5\_rc1 (released 28.12.18)
====================================

-  optimized GW-BSE code and integral engine
-  added closed shell DFT code which supports hybrid functionals
-  removed ctp dependency
-  atm no support for site energy calculation
-  CHELPG fit for ground and excited states
-  merged igwbse and idft into one calculator
-  added unit and integration test
-  improved geometry optimiser
-  replaced ublas with Eigen3
-  replaced boost serialisation with hdf5 files

Version 1.4.1 (released 02.09.17)
=================================

-  fix pkg-config file

Version 1.4 (released 29.10.16)
===============================

-  fixed a bug in gwbse
-  added missing copyright
-  cmake: fixed underlinking

Version 1.4\_rc1 (released 26.09.16)
====================================

-  include manual
-  an extension of the whole workflow from: electrons and holes, to
   singlet and triplet excitons
-  a fully functional GW-BSE code optimized for: molecular properties,
   including excited state geometry optimizsation
-  Inclusion of LIBXC to calculate Exchange correlation matrices
-  allowing interfacing GW-BSE with many quantum mechanical packages
-  support for ORCA DFT package
-  framework to use stochastic models to generate larger system
-  better performance of larger systems
-  new calculators: egwbse,igwbse,ewald,.....
-  support for intel mkl library and compilers for better performance
-  A periodic polarisation embedding: to calculate classical
   configuration energies without cutoffs
-  xtp\_update\_exciton to update state file to newest format
-  integration of moo and kmc into xtp for easier installation
-  kmc\_lifetime calculator to simulate exciton movement with lifetimes
-  partialcharges to extract atomic charges from qm calculation
-  renaming from ctp to xtp
-  many bugfixes

Version 1.3 (released XX.09.15)
===============================

-  new executables: ctp\_tools, ctp\_dump, ctp\_parallel,
   xtp\_testsuite, xtp\_update
-  ctp\_tools wraps light-weight tools that assist e.g. in generating
   the system mapping file
-  ctp\_dump extracts information from the state file to human-readable
   format
-  ctp\_parallel wraps heavy-duty job-based calculators: allows
   synchronization across processes
-  ctp\_testsuite provides an easy-to-use environment to run: selected
   tests, individual calculators
-  ctp\_update updates an existent state file to the current version
-  new calculators: edft, idft, pdb2map, xqmultipole, ...
-  edft / idft: provide interfaces to the GAUSSIAN, TURBOMOLE & NWCHEM
   package, using packages computes: couplings, internal energies,
   partial charges
-  pdb2map (generates a system mapping file from an input coordinate
   file)
-  xqmultipole computes classical configuration energies of: charged
   clusters embedded in a molecular environment
-  enhanced usability via the command-line help, tutorial & test-suite
-  a GUI tutorial assists with the first practical steps in using
   VOTCA-CTP
-  an extended and homogenized help system provides: short infos on
   individual calculator options from the command line

Version 1.0 (released 23.10.11)
===============================

-  parallel evaluation of site energies using: Thole model + GDMA -
   Tinker no longer required
-  much clearer input files (and many more checks for input errors)
-  most of calculators are parallel and can be used on a cluster
-  bug in zindo/ctp interface fixed
-  state file now contains: the atomistic trajectory, rigid fragments,
   conjugated segments
-  support for several MD frames
