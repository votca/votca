# How to prepare the files for a VOTCA XTP system

## Mapping file

Use the `-m` flag of `xtp_map` to get an initial mapping file.

split the initial molecule into more segments if useful. For small molecules this is rarely necessary. Polymers are more difficult.

For each segment you have to make a `.xyz` file.

## XYZ FILES

Create a `.xyz` file from the `xtp_map -m` output.

For molecules with long side chains it may be better to cut these off. Just remove the atoms from the `.xyz` file.

Choose either **ORCA**, **GAUSSIAN** or **NWCHEM** as your qmpackage.

Optimise this geometry for groundstate electron/hole/singlet/triplet depending on what you want to do later.
You can find input files for **ORCA** in the subfolders of [here](https://github.com/votca/xtp-tutorials/tree/electrostatics/LAMMPS/Thiophene/QC_FILES).
You can find input files for **GAUSSIAN** in [here](https://github.com/votca/xtp-tutorials/tree/electrostatics/GROMACS/Methane/QC_FILES).

Add the filepaths of the `.xyz` files to the mapping file.

## MPS FILES

Use the optimised geometries to calculate CHELPG charges and polarizabilities for groundstate electron/hole/singlet/triplet states.
You can find input files for **ORCA** in [here](https://github.com/votca/xtp-tutorials/tree/electrostatics/LAMMPS/Thiophene/MP_FILES).
You can find input files for **GAUSSIAN** in [here](https://github.com/votca/xtp-tutorials/tree/electrostatics/GROMACS/Methane/QC_FILES).

Use `xtp_tools -e log2mps` to create `.mps` files from the chelpg qmpackage output for each state.

Use `xtp_tools -e molpol` to create atomic polarizabilities for these `.mps` files.

Add these filepaths to the mapping file.

## REORGANISATION ENERGIES

Calculate groundstate energy in excited state geometries *UnX*.

Calculate excited state energies in groundstate geometries *UxN*.

You still have groudstate energy in groundstate geometry from the `.xyz` files *UnN*.

You still have excited state energy in  excited state geometry from the `.xyz` files *UxX*.

Add the differences between these energies to mapping file.

## FRAGMENTS

If your segment consists of multiple conjugated fragments ("rings") you might want to split it into multiple fragments, so in the mapping they can rotate with respect to each other.

For each fragment add the qmatoms and multipoles. If you have cut sidechains of replace these atoms by ":".

## LOCAL FRAME

Choose three atoms in each fragment. Normally just take the heaviest or those which are the most rigidly bond to each other.

## WEIGHTS

Unless you want to do something special just use atomic masses.

You are done now and can run the rest of VOTCA.





