For more detailed information about the changes see the history of the [repository](https://github.com/votca/xtp/commits/master).

## Version 1.5 _SuperVictor_ (released XX.01.19)
 * enable gitlab CI

## Version 1.5_rc3 (released 19.01.19)
 * travis: fixed bug in building tags

## Version 1.5_rc2 (released 16.01.19)
* fix parallel build of manual
* fix usage on inkscape on arm arch
* clean up namespace in header

## Version 1.5_rc1 (released 28.12.18)
* optimized GW-BSE code and integral engine
* added closed shell DFT code which supports hybrid functionals
* removed ctp dependency
* atm no support for site energy calculation
* CHELPG fit for ground and excited states
* merged igwbse and idft into one calculator
* added unit and integration test
* improved geometry optimiser
* replaced ublas with Eigen3
* replaced boost serialisation with hdf5 files

## Version 1.4.1 (released 02.09.17)
* fix pkg-config file

## Version 1.4 (released 29.10.16)
* fixed a bug in gwbse
* added missing copyright
* cmake: fixed underlinking 

## Version 1.4_rc1 (released 26.09.16)
* include manual
* an extension of the whole workflow from electrons and holes, to singelt and triplet excitons
* a fully functional GW-BSE code optimized for molecular properties, including excited state geometry optimizsation
* Inclusion of LIBXC to calculate Exchange correlation matrices, allowing interfacing GW-BSE with many quantum mechanical packages
* support for ORCA DFT package 
* framework to use stochastic models to generate larger system
* better performance of larger systems
* new calculators: egwbse,igwbse,ewald,.....
* support for intel mkl library and compilers for better performance
* Aperiodic polarisation embedding to calculate classical configuration energies without cutoffs
* xtp_update_exciton to update state file to newest format
* integration of moo and kmc into xtp for easier installation
* kmc_lifetime calculator to simulate exciton movement with lifetimes
* partialcharges to extract atomic charges from qm calculation
* renaming from ctp to xtp
* many bugfixes



## Version 1.3 (released XX.09.15)
* new executables: ctp_tools, ctp_dump, ctp_parallel, xtp_testsuite, xtp_update
* ctp_tools wraps light-weight tools that assist e.g. in generating the system mapping file
* ctp_dump extracts information from the state file to human-readable format
* ctp_parallel wraps heavy-duty job-based calculators, allowing for synchronization across processes
* ctp_testsuite provides an easy-to-use environment to run selected tests as well as individual calculators
* ctp_update updates an existent state file to the current version
* new calculators: edft, idft, pdb2map, xqmultipole, ...
* edft / idft (provide interfaces to the GAUSSIAN, TURBOMOLE & NWCHEM packages to compute couplings, internal energies and partial charges)
* pdb2map (generates a system mapping file from an input coordinate file)
* xqmultipole (computes classical configuration energies of charged clusters embedded in a molecular environment)
* enhanced usability via the command-line help, tutorial & test-suite
* a GUI tutorial assists with the first practical steps in using VOTCA-CTP
* an extended and homogenized help system provides short infos on individual calculator options from the command line

## Version 1.0 (released 23.10.11)
* parallel evaluation of site energies (Thole model + GDMA) - Tinker no longer required
* much clearer input files (and many more checks for input errors)
* most of calculators are parallel and can be used on a cluster
* bug in zindo/ctp interface fixed
* state file now contains the atomistic trajectory, rigid fragments and conjugates segments
* support for several MD frames
