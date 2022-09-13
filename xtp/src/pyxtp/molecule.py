"""Molecule representation.

API
---
.. autoclass:: Molecule

"""
import copy as cp
from pathlib import Path
from typing import List, Optional, Tuple, Union

import h5py
import numpy as np

from .utils import BOHR2ANG

Pathlike = Union[Path, str]


class Molecule:
    """Molecule definition."""

    def __init__(self):
        self.name = "molecule"
        self.elements = []
        self.coordinates = []
        self.gradient = None
        self.has_data = False
        self.has_xyz = False
        self.has_gradient = False
        self.hessian = None
        self.has_hessian = False

    def add_atom(self, element: str, x: float, y: float, z: float):
        """Add a single atom to the molecule."""
        self.elements.append(element)
        self.coordinates.append(np.array([x, y, z]))
        self.has_xyz = True

    def copy_and_displace(self, mol2, atomidx: int, coordidx: int, dr: float) -> None:
        """Create a new molecule displace by a dr factor."""
        if self.has_xyz:
            raise Exception("Molecule coordinates already defined!")

        # deep-copying elements and coordinates
        self.elements = cp.deepcopy(mol2.elements)
        self.coordinates = cp.deepcopy(mol2.coordinates)

        # displacing one atom in one direction as requested
        self.coordinates[atomidx][coordidx] += dr

    def print_xyz(self) -> None:
        """Print the molecule in xyz format."""
        for (element, coordinates) in zip(self.elements, self.coordinates):
            print(element, coordinates)

    def read_xyz_file(self, filename: Pathlike) -> None:
        """Read the molecular coordinates from a given file."""
        with open(filename, 'r') as handler:
            lines = handler.readlines()

        self.name = Path(filename).stem
        arr = [(row[0], np.array(row[1:], dtype=float)) for row in [
            x.split() for x in lines[2:]]]
        self.elements, self.coordinates = tuple(zip(*arr))
        self.has_xyz = True

    def write_xyz_file(self, filename: Pathlike):
        """Write the molecule in XYZ format."""
        atoms = "\n".join(f"{elem} {xyz[0]:.4f} {xyz[1]:.4f} {xyz[2]:.4f}" for elem, xyz in zip(
            self.elements, self.coordinates))
        mol = f"""{len(self.elements)}
{self.name} created by pyvotca writer
{atoms}
"""
        with open(filename, "w") as xyzfile:
            xyzfile.write(mol)

    def get_total_energy(self, kind: str, level: int, dynamic: bool = False) -> float:
        """Wrap call to individual total energy functions."""
        if kind == 'dft_tot':
            return self.get_dft_energy()
        elif kind == 'ks':
            return self.get_ks_total_energy(level)
        elif kind == 'qp_pert':
            return self.get_qp_total_energy(level)
        elif kind == 'qp_diag':
            return self.get_qp_total_energy(level)
        elif kind == 'bse_singlet' and not dynamic:
            return self.get_bse_singlet_total_energy(level)
        elif kind == 'bse_singlet' and dynamic:
            return self.get_bse_singlet_dynamic_total_energy(level)
        elif kind == 'bse_triplet' and not dynamic:
            return self.get_bse_triplet_total_energy(level)
        elif kind == 'bse_triplet' and dynamic:
            return self.get_bse_triplet_dynamic_total_energy(level)
        else:
            raise Exception(
                f'Energy of kind {kind} is not available!')

    def get_gradient(self):
        """Return the stored nuclear gradient in Hartree/Bohr."""
        if self.has_gradient:
            return self.gradient
        else:
            raise Exception(
                'Nuclear gradient not available!')

    def get_dft_energy(self):
        """Return the DFT total energy."""
        self.check_data()

        return self.DFTenergy

    def get_ks_total_energy(self, level=''):
        """Return the excited state KS total energy."""
        self.check_data()

        lumo = self.homo + 1

        total_energy = self.DFTenergy
        if (level < lumo):
            return(total_energy - self.ks_energies[level])
        elif level < len(self.ks_energies):
            return(total_energy + self.ks_energies[level])
        else:
            print("Requested KS level {} does not exist.")
            return 0.0

    def get_qp_total_energy(self, level=''):
        """Return the excited state QP total energy."""
        self.check_data()

        lumo = self.homo + 1

        total_energy = self.DFTenergy
        if (level < lumo):
            return(total_energy - self.qp_energies[level - self.qpmin])
        elif level < len(self.ks_energies):
            return(total_energy + self.qp_energies[level - self.qpmin])
        else:
            print("Requested QP level {} does not exist.")
            return 0.0

    def get_qp_diag_total_energy(self, level=''):
        """Return the excited state diag QP total energy."""
        self.check_data()

        lumo = self.homo + 1

        total_energy = self.DFTenergy
        if (level < lumo):
            return(total_energy - self.qp_energies_diag[level - self.qpmin])
        elif level < len(self.ks_energies):
            return(total_energy + self.qp_energies_diag[level - self.qpmin])
        else:
            print(f"Requested diag QP {level} does not exist.")
            return 0.0

    def get_bse_singlet_total_energy(self, level: int) -> float:
        """Return the excited state BSE Singlet total energy."""
        msg = f"Requested BSE singlet {level} does not exist."
        return self.check_and_read(level, "bse_singlet_energies", msg)

    def get_bse_triplet_total_energy(self, level: int) -> float:
        """Return the excited state BSE Singlet total energy."""
        msg = f"Requested BSE triplet {level} does not exist."
        return self.check_and_read(level, "bse_triplet_energies", msg)

    def get_bse_singlet_dynamic_total_energy(self, level: int) -> float:
        """Return the excited state BSE Singlet total energy."""
        msg = f"Requested dynamic BSE singlet {level} does not exist."
        return self.check_and_read(level, "bse_singlet_energies_dynamic", msg)

    def get_bse_triplet_dynamic_total_energy(self, level: int) -> float:
        """Return the excited state BSE Singlet total energy."""
        msg = f"Requested dynamic BSE triplet level {level} does not exist."
        return self.check_and_read(level, "bse_triplet_energies_dynamic", msg)

    def read_orb(self, orbfile: Pathlike) -> None:
        """Read data from the orb (HDF5) file."""
        with h5py.File(orbfile, 'r') as handler:
            orb = handler['QMdata']
            # get coordinates
            atoms = orb['qmmolecule']['qmatoms']
            # coordinates are stored in Bohr!
            arr = [(atom['element'][0].decode(), BOHR2ANG * np.array(
                [atom['posX'][0], atom['posY'][0], atom['posZ'][0]], dtype=float)) for atom in atoms]
            elements_in, coordinates_in = tuple(zip(*arr))

            if not self.has_xyz:
                self.elements = elements_in
                self.coordinates = coordinates_in
            else:
                self.check_molecule_integrity(elements_in, coordinates_in)

            self.has_xyz = True

            self.homo = int(orb.attrs['occupied_levels']) - 1
            self.DFTenergy = float(orb.attrs['qm_energy'])
            self.qp_energies = read_flatten_array(orb, 'QPpert_energies')
            self.qp_energies_diag, self.ks_energies, self.bse_singlet_energies, self.bse_triplet_energies = [
                read_flatten_array(orb, x, 'eigenvalues') for x in ('QPdiag', 'mos', 'BSE_singlet', 'BSE_triplet')]

            self.bse_singlet_energies_dynamic, self.bse_triplet_energies_dynamic = [
                read_flatten_array(orb, f"BSE_{x}_dynamic") for x in ("singlet", "triplet")]
            self.qpmin = int(orb.attrs['qpmin'])
            self.qpmax = int(orb.attrs['qpmax'])
            td = orb['transition_dipoles']
            self.transition_dipoles = np.array(
                [td[dset][()] for dset in td.keys()])
            self.has_data = True

    def check_molecule_integrity(self, other_elements: List[str], other_coordinates: List[np.ndarray]):
        """Compare the atoms from self with the one stored in the HDF5."""
        for k, (elem, coord, other_elem, other_coord) in enumerate(
                zip(self.elements, self.coordinates, other_elements, other_coordinates)):
            if elem != other_elem:
                raise Exception(
                    f'Element {elem} with index {k} in molecule differs from element {other_elem} in orb file!')

            if not np.allclose(coord, other_coord):
                raise Exception(
                    f'Molecular coordinates of element {k} {coord} differ from coordinates in orb file {other_coord}')

    def get_qp_corrections(self):
        self.check_data()

        qp_corrections = self.qp_energies -\
            self.ks_energies[self.qpmin:self.qpmin + len(self.qp_energies)]

        return qp_corrections.flatten()

    def get_oscillator_strengths(self, dynamic: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """Retrieve oscillator strenghts' values."""
        self.check_data()

        # get energies/oscillator strengths
        if dynamic:
            energy = self.bse_singlet_energies_dynamic
        else:
            energy = self.bse_singlet_energies
        osc = [(2. / 3.) * e * (t ** 2).sum()
               for e, t in zip(energy, self.transition_dipoles)]

        return energy, np.array(osc)

    def check_data(self):
        """Check that there is data in the molecule."""
        if not self.has_data:
            raise Exception("No energy has been stored!")

    def check_and_read(self, level: int, prop: str, msg: str) -> float:
        """Check that there is data available and retrieve it."""
        self.check_data()

        if level < len(getattr(self, prop)):
            return(self.DFTenergy + getattr(self, prop)[level])
        else:
            print(msg)
            return 0.0


def read_flatten_array(group: h5py.Group, key1: str, key2: Optional[str] = None):
    """Read an array from h5py handler and flatten it."""
    if key2 is None:
        arr = group[key1][()]
    else:
        arr = group[key1][key2][()]

    return arr.flatten()
