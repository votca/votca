import copy
from typing import List, Tuple, Optional, Union
from pyxtp.xtp import Pathlike
from pathlib import Path

import ase
import ase.io

class Molecule(ase.Atoms):
    """Molecule
    Represents a votca molecule. Inherits from ase.Atoms
    """

    AtomsLike = Union[ase.Atoms, List[ase.Atoms]]
    ListOfPositions = List[Tuple[float,float,float]]

    __slots__ = ("name", "elements", "coordinates", "gradient", "hessian")

    def __init__(self,
                 name: str = "molecule",
                 atoms: Optional[AtomsLike] = None,
                 elements: Optional[List[str]]=None,
                 coordinates: Optional[ListOfPositions]=None,
                 gradient=None,
                 hessian=None,
                 ):

        self.name = name

        if atoms is None:
            super().__init__(symbols=elements, positions=coordinates)
        else:
            super().__init__(atoms)

        self.elements = self.symbols
        self.coordinates = self.positions

        self.gradient = gradient
        self.hessian = hessian

    @property
    def has_data(self) -> bool:
        return False

    @property
    def has_coordinates(self) -> bool:
        return False

    @property
    def has_gradient(self) -> bool:
        return self.gradient is not None

    @property
    def has_hessian(self) -> bool:
        return self.hessian is not None

    def _add_ase_atom(self, atom: ase.Atom) -> None:
        self.append(atom)

    def add_atom(self, element: str, x: float, y: float, z: float):
        """Add a single atom to the molecule."""

        self._add_ase_atom(ase.Atom(element, (x,y,z)))

    def copy_and_displace(self,
                          mol2: "Molecule",
                          atomidx: int,
                          coordidx: int,
                          dr: float) -> "Molecule":
        """Create a new molecule displace by a dr factor."""

        new = copy.deepcopy(mol2)

        # displacing one atom in one direction as requested
        new.get_positions()[atomidx][coordidx] += dr

        return new

    @staticmethod
    def to_xyz_file(path_to_file: Pathlike, mol: "Molecule"):
        ase.io.write(path_to_file, images = mol, format = "xyz")

    @staticmethod
    def from_xyz_file(path_to_file: Pathlike):

        file_path = Path(path_to_file)
        if not file_path.exists():
            raise ValueError(f"{path_to_file} does not exist")

        atoms = ase.io.read(path_to_file, format="xyz")
        return Molecule(Path(path_to_file).stem, atoms=atoms)
