#! /usr/bin/env bash
# Remove all the output from the jupyter notebooks

files=("GROMACS/Methane/QMMM_GROMACS.ipynb" "LAMMPS/KMC_Thiophene/LAMMPS_KMC.ipynb" "LAMMPS/Thiophene/QMMM_LAMMPS.ipynb" "tools/dftgwbse_CH4/DFTGWBSE_ENERGY.ipynb" "tools/dftgwbse_CO_geoopt/DFTGWBSE_OPTIMIZATION.ipynb")

for x in ${files[@]}
do
    echo ${x}
    jupyter nbconvert --ClearOutputPreprocessor.enabled=True --inplace "${x}"
done

# Clean up temporal files
git clean -fdx
