import os
import numpy as np
from typing import Any, Dict, Union, List , Optional, Tuple
import h5py
from pathlib import Path
import copy as cp 

from ase.units import Hartree, Bohr
from ase import Atoms
from ase.calculators.calculator import Calculator, FileIOCalculator, Parameters, ReadError, equal
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
        
        self.has_data = False
        self.has_gradient = False
        self.jobdir = './'
        self.hdf5_filename = ''
        self.logfile = ''
        # self.set_from_options(options)
    
    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2, ...).

        A dictionary containing the parameters that have been changed
        is returned.

        Subclasses must implement a set() method that will look at the
        chaneged parameters and decide if a call to reset() is needed.
        If the changed parameters are harmless, like a change in
        verbosity, then there is no need to call reset().
        """

        translation_dic = {'xc': 'dftpackage/functional',
                           'basis': 'dftpackage/basisset'}

        changed_parameters = {}

        for key, value in kwargs.items():
            
            if key in translation_dic:
                key = translation_dic[key]
                
            oldvalue = self.parameters.get(key)
            if key not in self.parameters or not equal(value, oldvalue):
                changed_parameters[key] = value
                self.parameters[key] = value
                
                split_key = key.split('/')
                element = self.options
                for k in split_key[:-1]:
                    element = element.__getattribute__(k)
                element.__setattr__(split_key[-1], value)
                
        if self.discard_results_on_any_change and changed_parameters:
            self.reset()
            
        return changed_parameters    

    def set_from_options(self, options: XTPOptions):
        """Set the options from an XTPOptions instance"""
        self.options = options
        if options is not None:
            opt_dict = options.__todict__()
            changed_parameters = self.set(opt_dict)
            if changed_parameters:
                self.reset()
                
    def set_atoms(self, atoms):
        """set atomic positions"""
        self.atoms = atoms
                
    def calculate(self, atoms=None, name=None):
        """Calculate things."""
        
        Calculator.calculate(self, atoms)
        atoms = self.atoms

        # write input files
        if name is None:
            xyzname = self.atoms.get_chemical_formula()
        else:
            xyzname = name
            
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
            self.hdf5_filename = f'{self.jobdir}{self.jobname}.orb'
            os.replace(f'{xyzname}.orb', self.hdf5_filename)
        else:
            self.hdf5_filename = f'{xyzname}.orb'
            self.logfile = 'dftgwbse.log'
            
        # Reads energies from an existing HDF5 orb file
        self.read_hdf5_data(self.hdf5_filename)
   
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
        
    def read_hdf5_data(self, hdf5_filename: Pathlike) -> None:
        """Read data from the orb (HDF5) file."""
        with h5py.File(hdf5_filename, 'r') as handler:
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

    def read_forces_from_logfile(self, logfile: Pathlike) -> None:
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

    def calculate_displaced_geometries(self, eps = 0.001):
        """Run a VOTCA simulation for every displacement of the atomic positions."""
          
        def copy_and_displace_atoms(atoms: Atoms, idx_atom: int,  
                                    idx_coordinate: int, eps: float):
            """Make a copy of the atoms and dispace atoim idx_atom in direction idx_coordinate by eps

            Args:
                atoms (Atoms): atoms object
                idx_atom (int): index of thew atom to  move
                idx_coordinate (int): index of the coordinate to move
                eps (float): displacement
            """
            atoms_displaced = atoms.copy()
            positions = atoms_displaced.get_positions()
            positions[idx_atom][idx_coordinate] += eps
            atoms_displaced.set_positions(positions)
            return atoms_displaced           
            
        # how many atoms
        natoms = len(self.atoms)
        directions = [-1.0, 1.0]
        self.atoms_displaced = dict()
        self.atoms_eps = eps 
        
        for atom in range(natoms):
            for coordinate in range(3):
                for direction in directions:
                    
                    # get a unique name
                    name = self._generate_gradient_name(self.atoms.get_chemical_formula(),
                                                        atom,
                                                        direction,
                                                        coordinate)
                    
                    # get displaced molecule
                    self.atoms_displaced[name] = copy_and_displace_atoms(self.atoms,
                                                              atom,
                                                              coordinate,
                                                              float(direction)* eps * BOHR2ANG)
                    
                    # make a new xtp wrapper for this one
                    calc = xtp()
                    calc.set_from_options(self.options)
                    self.atoms_displaced[name].calc = calc
                    
                    # run this
                    self.atoms_displaced[name].calc.calculate(name=name)
                    
    def calculate_numerical_gradient(self, kind: str, energy_level: int) -> np.ndarray:
        """Computes the gradient for a particle/excitation kind expecting
        all displaced calculation to be available.

        Parameters
        ----------
        kind of particle/excitation (choices: BSE_singlet, BSE_triplet, QPdiag, QPpert and dft_tot)
          and optionally the energy level if not provided all energy levels will be returned

        """

        # how many atoms
        natoms = len(self.atoms)
        
        # store gradient in xtp.mol object
        self.results['forces'] = np.zeros((natoms, 3))

        directions = [-1.0, 1.0]
        for atom in range(natoms):
            for coordinate in range(3):
                energy_plus = 0.0
                energy_minus = 0.0
                for direction in directions:
                    # get energy for displaced molecules
                    name = self._generate_gradient_name(self.atoms.get_chemical_formula(),
                                                        atom,
                                                        direction,
                                                        coordinate)
                    atoms_displaced = self.atoms_displaced[name]
                    
                    if direction > 0:
                        energy_plus = atoms_displaced.calc.get_total_energy(
                            kind, energy_level)
                    else:
                        energy_minus = atoms_displaced.calc.get_total_energy(
                            kind, energy_level)

                self.results['forces'][atom, coordinate] = (
                    energy_plus - energy_minus) / (2.0 * self.atoms_eps)

        self.has_gradient = True

    def get_gradient(self, kind: str, energy_level: int) -> np.ndarray:
        """Retrieve the gradient."""
        return self.calculate_numerical_gradient(kind, energy_level)

    @staticmethod
    def _generate_gradient_name(self, name: str, atom: int, dir: float, coord: int) -> str:
        """Generate a name for the gradient calculation."""
        return f"{name}_{atom}_{dir}_{coord}"

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


