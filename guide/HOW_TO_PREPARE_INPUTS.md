Mapping file

use the -m flag of xtp_map to get an initial mapping file

split the initial molecule into more segments if useful. For small molecules this is rarely necessary. Polymers are more difficult.

for each segment

XYZ FILES

create a .xyz file from the xtp_map -m output

For molecules with long side chains it may be better to cut these off. Just remove the atoms from the .xyz file

choose either ORCA, GAUSSIAN or NWCHEM as your qmpackage


optimise this geometry for groundstate electron/hole/singlet/triplet depending on what you want to do later
you can find input files for ORCA in the thiophene QC_FILES folder in Thiophene
you can find input files for gaussian in the methane QC_FILES folder in Methane

add these filepaths to the mapping file

MPS FILES

use the optimised geometries to calculate chelpg charges and polarizabilities for groundstate electron/hole/singlet/triplet states
you can find input files for ORCA in the thiophene MP_FILES folder in Thiophene
you can find input files for gaussian in the methane MP_FILES folder in Methane

use xtp_tools -e log2mps to create .mps files from the chelpg qmpackage output for each state

use xtp_tools -e molpol to create atompolarizabilities for these .mps files

add these filepaths to the mapping file

REORGANISATION ENERGIES

calculate groundstate energy in excited state geometries UnX

calculate excited state energies in groundstate geometries UxN

you still have groudstate energy in groundstate geometry from the XYZ files UnN

you still have excited state energy in  excited state geometry from the XYZ files UxX

add the differences between these energies to mapping file

FRAGMENTS

if your segment consists of multiple conjugated fragments ("rings") you might want to split it into multiple fragments, so in the mapping they can rotate with respect to each other

for each fragment add the qmatoms and multipoles, if you have cut sidechains of replace these atoms by ":"

LOCAL FRAME

Choose three atoms in each fragment. Normally just take the heaviest or those which are the most rigidly bond to each other.

WEIGHTS

Unless you want to do something special just use atomic masses

you are done now and can run the rest of VOTCA





