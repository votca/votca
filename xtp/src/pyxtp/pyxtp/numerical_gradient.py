"""Numerical gradient.

API
---
.. autoclass:: NumericalGradient

"""

import numpy as np

from .molecule import Molecule
from .utils import BOHR2ANG
from .dftgwbse import DFTGWBSE

__all__ = ["NumericalGradient"]


class NumericalGradient:
    """Compue the gradient numerically using XTP."""

    def __init__(self, xtp: DFTGWBSE, dr: float = 0.001, path_to_simulations: str = './gradient/'):
        self.xtp = xtp
        self.dr = dr
        self.path = path_to_simulations

    def gen_name(self, name: str, atom: int, dir: float, coord: int) -> str:
        """Generate a name for the gradient calculation."""
        return f"{name}_{atom}_{dir}_{coord}"

    def run_permut(self):
        """Run a VOTCA simulation for every displacement of the electric field with strength dE."""
        # Initial structure expected in the mol object of xtp

        # how many atoms
        natoms = len(self.xtp.mol.elements)

        directions = [-1.0, 1.0]
        for atom in range(natoms):
            for coordinate in range(3):
                for direction in directions:
                    # get displaced molecule
                    mol_displaced = Molecule()
                    mol_displaced.copy_and_displace(
                        self.xtp.mol, atom, coordinate, float(direction) * self.dr * BOHR2ANG)
                    name = self.gen_name(
                        mol_displaced.name, atom, direction, coordinate)
                    # make a new xtp wrapper for this one
                    xtp_displaced = DFTGWBSE(mol_displaced, threads=self.xtp.threads,
                                             options=self.xtp.options, jobname=name, jobdir=self.path)
                    # run this
                    xtp_displaced.run()

    def calc_gradient(self, kind: str, energy_level: int) -> np.ndarray:
        """Computes the gradient for a particle/excitation kind expecting
        all displaced calculation to be available.

        Parameters
        ----------
        kind of particle/excitation (choices: BSE_singlet, BSE_triplet, QPdiag, QPpert and dft_tot)
          and optionally the energy level if not provided all energy levels will be returned

        Returns
        -------
        Numpy array of nuclear gradient stored in molecule object.     

        """

        # how many atoms
        natoms = len(self.xtp.mol.elements)
        # store gradient in xtp.mol object
        self.xtp.mol.gradient = np.zeros((natoms, 3))

        directions = [-1.0, 1.0]
        for atom in range(natoms):
            for coordinate in range(3):
                energy_plus = 0.0
                energy_minus = 0.0
                for direction in directions:
                    # get energy for displaced molecules
                    mol_displaced = Molecule()
                    name = self.gen_name(
                        mol_displaced.name, atom, direction, coordinate)
                    orbname = self.path + name + '.orb'
                    mol_displaced.read_orb(orbname)
                    if direction > 0:
                        energy_plus = mol_displaced.get_total_energy(
                            kind, energy_level)
                    else:
                        energy_minus = mol_displaced.get_total_energy(
                            kind, energy_level)

                self.xtp.mol.gradient[atom, coordinate] = (
                    energy_plus - energy_minus) / (2.0 * self.dr)

        self.xtp.mol.has_gradient = True
        return self.xtp.mol.gradient

    def get_gradient(self, kind: str, energy_level: int) -> np.ndarray:
        """Retrieve the gradient."""
        return self.calc_gradient(kind, energy_level)
