For more detailed information about the changes see the history of the [repository](https://github.com/votca/xtp/commits/master).

## Version 1.3 (released XX.09.15)
* new executables: xtp_tools, xtp_dump, xtp_parallel, xtp_testsuite, xtp_update
* xtp_tools wraps light-weight tools that assist e.g. in generating the system mapping file
* xtp_dump extracts information from the state file to human-readable format
* xtp_parallel wraps heavy-duty job-based calculators, allowing for synchronization across processes
* xtp_testsuite provides an easy-to-use environment to run selected tests as well as individual calculators
* xtp_update updates an existent state file to the current version
* new calculators: edft, idft, pdb2map, xqmultipole, ...
* edft / idft (provide interfaces to the GAUSSIAN, TURBOMOLE & NWCHEM packages to compute couplings, internal energies and partial charges)
* pdb2map (generates a system mapping file from an input coordinate file)
* xqmultipole (computes classical configuration energies of charged clusters embedded in a molecular environment)
* enhanced usability via the command-line help, tutorial & test-suite
* a GUI tutorial assists with the first practical steps in using VOTCA-XTP
* an extended and homogenized help system provides short infos on individual calculator options from the command line

## Version 1.0 (released 23.10.11)
* parallel evaluation of site energies (Thole model + GDMA) - Tinker no longer required
* much clearer input files (and many more checks for input errors)
* most of calculators are parallel and can be used on a cluster
* bug in zindo/xtp interface fixed
* state file now contains the atomistic trajectory, rigid fragments and conjugates segments
* support for several MD frames
