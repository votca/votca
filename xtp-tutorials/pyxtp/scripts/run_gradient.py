#!/usr/bin/env python
"""Example to perform a gradient calculation."""
from pyxtp import xtp
from ase import Atoms
import numpy as np

def run_gradient() -> np.ndarray:
    """Compute gradient."""
    # define a molecule
    atoms = Atoms('CO', positions=([0,0,0],[1.4,0,0]))

    # define the calculator
    calc = xtp(nthreads=2)
    calc.select_force(energy='singlets', level=0, dynamic=False)
    
    # this allows to change all options
    # calc.options.dftpackage.functional = 'PBE'
    calc.options.dftpackage.basisset = 'def2-svp'
    calc.options.dftpackage.auxbasisset = 'aux-def2-svp'

    # attach the calculator
    atoms.calc = calc

    # run for the molecule in its geometry
    return atoms.get_forces()

if __name__ == "__main__":
    run_gradient()
