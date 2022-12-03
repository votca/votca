import os
import numpy as np
from typing import Any, Dict, Union, List , Optional, Tuple
import h5py
from pathlib import Path

from ase.units import Hartree, Bohr
from ase import Atoms
from ase.calculators.calculator import Calculator, FileIOCalculator, Parameters, ReadError
from ..capture_standard_output import capture_standard_output
from pyxtp import xtp_binds
from ..options import XTPOptions
from ..molecule import Molecule
from ..dftgwbse import DFTGWBSE 
from ..utils import BOHR2ANG

Pathlike = Union[Path, str]

def load_default_parameters():
    """Create a dictionary of default parameters"""
    options = XTPOptions()
    options._fill_default_values()
    return options.__todict__()

class xtp(Calculator):
    """This the ASE-calculator for doing XTP calculation."""

    implemented_properties = ['energy', 'forces', 'singlets',
                              'triplets', 'qp', 'ks', 
                              'qp_pert', 'transition_dipoles']

    default_parameters: Dict[str, Any] = load_default_parameters()

    def __init__(self, 
                 restart=None, 
                 *,
                 label=None, 
                 atoms=None,
                 nthreads=1,
                 **kwargs):
        
        Calculator.__init__(self, restart=restart, label=label, atoms=atoms, **kwargs)
        self.options = XTPOptions()
        self.nthreads = nthreads
        
        self.hdf5_data = None
        self.has_data = False
        self.has_data = False
        self.has_gradient = False
        self.jobdir = './'
        self.hdf5_data = ''
        self.logfile = ''
        # self.set_from_options(options)
            

    def set_from_options(self, options: XTPOptions):
        """Set the options from an XTPOptions instance"""
        if options is not None:
            opt_dict = options.__todict__()
            changed_parameters = self.set(opt_dict)
            if changed_parameters:
                self.reset()
                
    def set_atoms(self, atoms):
        """set atomic positions"""
        self.atoms = atoms
                
    def calculate(self, atoms=None):
        """Calculate things."""
        
        Calculator.calculate(self, atoms)
        atoms = self.atoms

        # write input files
        xyzname = 'molecule'
        xyz_filename = xyzname + '.xyz'
        input_filename = 'dftgwbse.xml'
        self.atoms.write(xyz_filename)
        self.options.job_name = xyzname
        self.options._write_xml(input_filename)

        """ Runs VOTCA and moves results a job folder, if requested """
        if not Path(self.jobdir).exists():
            os.makedirs(self.jobdir)

        # current dftgwbse file
        path_dftgwbse = (Path(input_filename)).absolute().as_posix()
        
        # Call the tool and capture the standard output
        output = capture_standard_output(
            xtp_binds.call_tool, "dftgwbse", self.nthreads, path_dftgwbse)
        
        with open(xyzname + ".out", "w") as handler:
            handler.write(output)

        # copy orbfile, if jobdir is not default
        if (self.jobdir != "./"):
            self.hdf5_data = f'{self.jobdir}{self.jobname}.orb'
            os.replace(f'{xyzname}.orb', self.hdf5_data)
        else:
            self.hdf5_data = f'{xyzname}.orb'
            self.logfile = 'dftgwbse.log'
            
        # Reads energies from an existing HDF5 orb file
        self.read_hdf5_data(self.hdf5_data)
   
        # create the results        
        self.results = {
            'energy': self.DFTenergy,
            'singlets': self.bse_singlet_energies,
            'triplets': self.bse_triplet_energies,
            'ks': self.ks_energies,
            'qp': self.qp_energies_diag,
            'transition_dipoles': self.transition_dipoles,
            'qp_pert': self.qp_energies,
            # 'forces': self.atomic_forces * Hartree / Bohr   
        }
        
    def read_hdf5_data(self, orbfile: Pathlike) -> None:
        """Read data from the orb (HDF5) file."""
        with h5py.File(orbfile, 'r') as handler:
            orb = handler['QMdata']
            # get coordinates
            atoms = orb['qmmolecule']['qmatoms']
            # coordinates are stored in Bohr!
            arr = [(atom['element'][0].decode(), BOHR2ANG * np.array(
                [atom['posX'][0], atom['posY'][0], atom['posZ'][0]], dtype=float)) for atom in atoms]
            elements_in, coordinates_in = tuple(zip(*arr))

            if self.atoms is None:
                self.atoms = Atoms(elements_in, positions=coordinates_in)
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

    def read_forces(self, logfile: Pathlike) -> None:
        """Read Forces from VOTCA logfile."""
        fil = open(logfile, 'r', encoding='utf-8')
        lines = fil.readlines()
        fil.close()
        getgrad = "no"
        for i, line in enumerate(lines):
            if line.find('ENGRAD') >= 0:
                getgrad = "yes"
                gradients = []
                tempgrad = []
                continue
            if getgrad == "yes" and "Saving" not in line:
                if "Total" not in line:
                    grad = line.split()
                    tempgrad.append(float(grad[3]))
                    tempgrad.append(float(grad[4]))
                    tempgrad.append(float(grad[5]))
                    gradients.append(tempgrad)
                    tempgrad = []
                else:
                    energy = float(line.split()[-2])
            if 'Saving' in line:
                getgrad = "no"
        self.atomic_forces = -np.array(gradients) #* Hartree / Bohr

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

    def check_molecule_integrity(self, other_elements: List[str], other_coordinates: List[np.ndarray]):
        """Compare the atoms from self with the one stored in the HDF5."""
        for k, (elem, coord, other_elem, other_coord) in enumerate(
                zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions(), 
                    other_elements, other_coordinates)):
            if elem != other_elem:
                raise Exception(
                    f'Element {elem} with index {k} in molecule differs from element {other_elem} in orb file!')

            if not np.allclose(coord, other_coord):
                raise Exception(
                    f'Molecular coordinates of element {k} {coord} differ from coordinates in orb file {other_coord}')


def read_flatten_array(group: h5py.Group, key1: str, key2: Optional[str] = None):
    """Read an array from h5py handler and flatten it."""
    if key2 is None:
        arr = group[key1][()]
    else:
        arr = group[key1][key2][()]

    return arr.flatten()


