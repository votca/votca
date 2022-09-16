#!/usr/bin/env python
"""Example to compute energies using XTP."""
from pyxtp import DFTGWBSE, Molecule, Visualization

def run_energy(save_figure: bool = False):
    """Run energy workflow."""
    # define a molecule
    mol = Molecule()

    # make it by hand
    mol.add_atom("C", 0, 0, 0)
    mol.add_atom("O", 1.2, 0, 0)

    # or read it from existing file
    # mol.readXYZfile('CO.xyz')

    # get a DFTGWBSE object
    dft = DFTGWBSE(mol)

    # change basis sets to a smaller one
    dft.options.dftpackage.basisset = 'def2-svp'
    dft.options.dftpackage.auxbasisset = 'aux-def2-svp'

    # run for the molecule
    dft.run()

    # only needed, if no run was performed but an existing HDF5 is read
    dft.mol.read_orb('pyvotca/examples/example.orb')

    # Getting the plotting functions
    viz = Visualization(dft.mol, save_figure=save_figure)
    # plotting QP corrections
    # viz.plot_qp_corrections()
    # plotting absorption spectrum
    viz.plot_absorption_gaussian()


if __name__ == "__main__":
    run_energy()
