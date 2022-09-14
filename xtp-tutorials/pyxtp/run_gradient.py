#!/usr/bin/env python
"""Example to perform a gradient calculation."""
import numpy as np

from pyxtp import DFTGWBSE, Molecule, NumericalGradient


def run_gradient() -> np.ndarray:
    """Compute gradient."""
    # define a molecule
    mol = Molecule()

    # make it by hand
    mol.add_atom("C", 0, 0, 0)
    mol.add_atom("O", 1.3, 0.0, 0.0)

    # get a DFTGWBSE object
    votca = DFTGWBSE(mol)
    # this allows to change all options
    # votca.options['functional'] = 'PBE'
    votca.options.basisset = 'def2-svp'
    votca.options.auxbasisset = 'aux-def2-svp'

    # run for the molecule in its geometry
    votca.run()

    # calculate a DFT-GWBSE gradient at the geometry
    grad = NumericalGradient(votca, dr=0.001)
    grad.run_permut()
    grad.get_gradient('BSE_singlet', 0)

    print(mol.get_gradient())
    return mol.get_gradient()


if __name__ == "__main":
    run_gradient()
