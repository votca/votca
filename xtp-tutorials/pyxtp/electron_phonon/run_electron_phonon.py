#!/usr/bin/env python
"""Phonon electron example."""
from pathlib import Path
from typing import Tuple

import numpy as np

from pyvotca import Electronphonon, Molecule, Orca

PATH_EXAMPLES = Path("pyvotca/examples/electron_phonon")


def run_electron_phonon(show: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Run an electron phonon simulation."""
    # define a molecule
    mol = Molecule()

    # load xyz
    mol.read_xyz_file(PATH_EXAMPLES / 'NPB.xyz')
    orca = Orca(mol)

    # read gradient from orca
    orca.read_gradient(PATH_EXAMPLES / 'NPB+.engrad')

    # read Hessian from orca
    orca.read_hessian(PATH_EXAMPLES / 'NPB.hess')

    # calculate el-ph couplings and plot
    ep = Electronphonon()
    return ep.calculate_electron_phonon_couplings(mol.elements, mol.hessian, mol.gradient, plot=show)


if __name__ == "__main__":
    run_electron_phonon()
