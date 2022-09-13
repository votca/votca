"""Orca reader/writer module."""
from pathlib import Path
from typing import Union

from ..molecule import Molecule
from ..parsers.orca_parsers import parse_gradient, parse_hessian
from ..utils import BOHR2ANG

PathLike = Union[str, Path]

__all__ = ["Orca"]


class Orca:
    def __init__(self, mol: Molecule):
        self.mol = mol

    def read_gradient(self, gradient_file: PathLike):
        """Read the nuclear gradient from an orca engrad file."""
        # read the gradient from an ORCA engrad calculation
        self.mol.gradient = parse_gradient(gradient_file)
        self.mol.has_gradient = True

    def read_hessian(self, hessian_file: PathLike):
        """Read Hessian from the orca .hess file."""
        self.mol.hessian = parse_hessian(hessian_file)
        self.mol.has_hessian = True

    def write_gradient(self, energy: float):
        """Write gradient in Orca format."""
        s = f"{energy:.18g}\n"
        for gradient in self.mol.gradient:
            gs = [g / BOHR2ANG for g in gradient]
            s += f" {gs[0]:.18g} {gs[1]:.18g} {gs[2]:.18g}\n"

        with open("system.extcomp.out", "w") as orcaout:
            orcaout.write(s)
