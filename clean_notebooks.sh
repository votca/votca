#! /bin/bash
# Remove all the output from the jupyter notebooks

files=("GROMACS/Methane/QMMM_GROMACS.ipynb" "LAMMPS/KMC_Thiophene/LAMMPS_KMC.ipynb" "LAMMPS/Thiophene/QMMM_LAMMPS.ipynb" "tools/dftgwbse_CH4/DFTGWBSE_ENERGY.ipynb" "tools/dftgwbse_CO_geoopt/DFTGWBSE_OPTIMIZATION.ipynb")
dirs=("GROMACS/Methane" "LAMMPS/KMC_Thiophene" "LAMMPS/Thiophene" "tools/dftgwbse_CH4" "tools/dftgwbse_CO_geoopt")

for x in "${files[*]}"
do
    jupyter nbconvert --ClearOutputPreprocessor.enabled=True --inplace $x
done

for d in "${dirs[*]}"
do
    git clean -f
done
