#! /usr/bin/env bash
# Remove all the output from the jupyter notebooks

files=("thiophene/MAPPING.ipynb" "thiophene/CHARGE_TRANSPORT.ipynb" "thiophene/EXCITATIONS.ipynb" "tools/dftgwbse_CH4/DFTGWBSE_ENERGY.ipynb" "tools/transport_lattice/TRANSPORT.ipynb" "tools/stochastic_network/STOCHASTIC_NETWORK.ipynb")

for x in ${files[@]}
do
    echo ${x}
    jupyter nbconvert --ClearOutputPreprocessor.enabled=True --inplace "${x}"
done

# Clean up temporal files
git clean -fdx
