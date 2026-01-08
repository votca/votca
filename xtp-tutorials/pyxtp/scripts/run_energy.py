#!/usr/bin/env python3
"""Example to compute energies using XTP."""
from pyxtp import xtp, Visualization
from ase import Atoms

def run_energy(save_figure: bool = False):
    """Run energy workflow."""
    # define a molecule
    atoms = Atoms('CO', positions=([0,0,0],[1.4,0,0]))

    # define the calculator
    calc = xtp(nthreads=2)

    # change basis sets to a smaller one
    calc.options.dftpackage.basisset = 'def2-svp'
    calc.options.dftpackage.auxbasisset = 'aux-def2-svp'
    calc.options.logging_file = 'CO_energy.log'

    # attach the calculator
    atoms.calc = calc

    # run the calculations
    atoms.get_potential_energy()

    # only needed, if no run was performed but an existing HDF5 is read
    # atoms.read('pyvotca/examples/example.orb')

    # Getting the plotting functions
    viz = Visualization(atoms, save_figure=save_figure)
    # plotting QP corrections
    # viz.plot_qp_corrections()
    # plotting absorption spectrum
    viz.plot_absorption_gaussian()


if __name__ == "__main__":
    run_energy()
