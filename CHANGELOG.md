For more detailed information about the changes see the history of the
[repository](https://github.com/votca/xtp/commits/master).

## Version 1.4.1 (released 02.09.17)

* fix pkg-config file

## Version 1.4 (released 29.10.16)

* fixed a bug in gwbse
* added missing copyright
* cmake: fixed underlinking

## Version 1.4_rc1 (released 26.09.16)

* include manual
* an extension of the whole workflow from:
 * electrons and holes
 * to singlet and triplet excitons
* a fully functional GW-BSE code optimized for:
 * molecular properties, including excited state geometry optimizsation
* Inclusion of LIBXC to calculate Exchange correlation matrices
 * allowing interfacing GW-BSE with many quantum mechanical packages
* support for ORCA DFT package
* framework to use stochastic models to generate larger system
* better performance of larger systems
* new calculators: egwbse,igwbse,ewald,.....
* support for intel mkl library and compilers for better performance
* A periodic polarisation embedding:
 * to calculate classical configuration energies without cutoffs
* xtp_update_exciton to update state file to newest format
* integration of moo and kmc into xtp for easier installation
* kmc_lifetime calculator to simulate exciton movement with lifetimes
* partialcharges to extract atomic charges from qm calculation
* renaming from ctp to xtp
* many bugfixes

## Version 1.3 (released XX.09.15)

* new executables: ctp_tools, ctp_dump, ctp_parallel, xtp_testsuite, xtp_update
* ctp_tools wraps light-weight tools that assist e.g.
 * in generating the system mapping file
* ctp_dump extracts information from the state file to human-readable format
* ctp_parallel wraps heavy-duty job-based calculators:
 * allows synchronization across processes
* ctp_testsuite provides an easy-to-use environment to run:
 * selected tests
 * individual calculators
* ctp_update updates an existent state file to the current version
* new calculators: edft, idft, pdb2map, xqmultipole, ...
* edft / idft:
 * provide interfaces to the GAUSSIAN, TURBOMOLE & NWCHEM packages
 * using packages computes: 
  * couplings
  * internal energies
  * partial charges
* pdb2map (generates a system mapping file from an input coordinate file)
* xqmultipole computes classical configuration energies of:
 * charged clusters embedded in a molecular environment)
* enhanced usability via the command-line help, tutorial & test-suite
* a GUI tutorial assists with the first practical steps in using VOTCA-CTP
* an extended and homogenized help system provides:
 * short infos on individual calculator options from the command line

## Version 1.0 (released 23.10.11)

* parallel evaluation of site energies using:
 * Thole model + GDMA - Tinker no longer required
* much clearer input files (and many more checks for input errors)
* most of calculators are parallel and can be used on a cluster
* bug in zindo/ctp interface fixed
* state file now contains:
 * the atomistic trajectory
 * rigid fragments
 * conjugated segments
* support for several MD frames
