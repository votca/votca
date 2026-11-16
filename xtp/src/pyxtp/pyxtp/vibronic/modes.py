"""Normal modes from a Cartesian Hessian.

Units follow ASE: energies in eV, lengths in Angstrom, masses in amu,
Hessians in eV/Angstrom^2 and gradients in eV/Angstrom.  Mass-weighted
quantities carry the corresponding powers of amu^(1/2).

Translations and rotations are projected out exactly (Eckart conditions
at the reference geometry): the Hessian is diagonalized in the orthogonal
complement of the 6 (5 for linear molecules) rigid-body directions, so
there are always 3N-6 (3N-5) vibrational modes and no frequency threshold
is needed to remove them.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
from ase import Atoms, units

Pathlike = Union[str, Path]

#: hbar*omega in eV for an eigenvalue of 1 eV/(Angstrom^2 amu)
HBAR_OMEGA_EV = units._hbar * 1e10 / np.sqrt(units._e * units._amu)
#: 1 eV in cm^-1
EV2INVCM = units._e / (units._hplanck * units._c * 100.0)
#: Boltzmann constant in eV/K
KB_EV = units.kB


def rigid_body_vectors(positions: np.ndarray, masses: np.ndarray,
                       tol: float = 1e-6) -> np.ndarray:
    """Orthonormal mass-weighted translation/rotation vectors (3N x n_rb).

    n_rb is 6, or 5 for a linear molecule (3 for a single atom).
    """
    natoms = len(masses)
    sqm = np.sqrt(masses)
    com = (masses[:, None] * positions).sum(axis=0) / masses.sum()
    r = positions - com
    vecs = []
    for k in range(3):
        t = np.zeros((natoms, 3))
        t[:, k] = sqm
        vecs.append(t.ravel())
    for k in range(3):
        axis = np.zeros(3)
        axis[k] = 1.0
        rot = np.cross(axis, r) * sqm[:, None]
        vecs.append(rot.ravel())
    d = np.array(vecs).T
    # orthonormalize and drop dependent directions (linear molecules)
    u, s, _ = np.linalg.svd(d, full_matrices=False)
    return u[:, s > tol * s.max()]


@dataclasses.dataclass
class NormalModes:
    """Harmonic vibrations of a molecule at a reference geometry.

    Attributes:
        symbols: chemical symbols (N)
        positions: reference geometry (N, 3), Angstrom
        masses: atomic masses (N), amu
        hessian: Cartesian Hessian (3N, 3N), eV/Angstrom^2
        eigenvalues: mass-weighted force constants of the vibrations
            (n_vib), eV/(Angstrom^2 amu); negative for imaginary modes
        modes: orthonormal mass-weighted eigenvectors (3N, n_vib)
    """

    symbols: Sequence[str]
    positions: np.ndarray
    masses: np.ndarray
    hessian: np.ndarray
    eigenvalues: np.ndarray
    modes: np.ndarray

    # ------------------------------------------------------------------
    @classmethod
    def from_hessian(cls, atoms: Atoms, hessian: np.ndarray,
                     masses: Optional[np.ndarray] = None) -> "NormalModes":
        """Diagonalize the projected mass-weighted Hessian.

        Args:
            atoms: reference geometry (masses from atoms.get_masses()
                unless given)
            hessian: Cartesian Hessian (3N, 3N) in eV/Angstrom^2; it is
                symmetrized
            masses: optional masses (amu), e.g. for isotopes
        """
        h = np.asarray(hessian, dtype=float)
        n3 = 3 * len(atoms)
        if h.shape != (n3, n3):
            raise ValueError(f"Hessian has shape {h.shape}, expected "
                             f"({n3}, {n3})")
        h = 0.5 * (h + h.T)
        m = np.asarray(atoms.get_masses() if masses is None else masses,
                       dtype=float)
        pos = atoms.get_positions()
        inv_sqm = np.repeat(1.0 / np.sqrt(m), 3)
        hm = h * inv_sqm[:, None] * inv_sqm[None, :]
        # basis of the vibrational subspace = complement of rigid motions
        d = rigid_body_vectors(pos, m)
        q, _ = np.linalg.qr(np.hstack([d, np.eye(n3)]))
        b = q[:, d.shape[1]:n3]
        # remove the rigid-body components exactly
        b -= d @ (d.T @ b)
        b, _ = np.linalg.qr(b)
        hv = b.T @ hm @ b
        w, u = np.linalg.eigh(0.5 * (hv + hv.T))
        return cls(symbols=list(atoms.get_chemical_symbols()),
                   positions=pos.copy(), masses=m.copy(), hessian=h,
                   eigenvalues=w, modes=b @ u)

    # ------------------------------------------------------------------
    @property
    def natoms(self) -> int:
        return len(self.masses)

    @property
    def nmodes(self) -> int:
        return len(self.eigenvalues)

    @property
    def energies(self) -> np.ndarray:
        """Vibrational quanta hbar*omega in eV (negative: imaginary)."""
        return np.sign(self.eigenvalues) * HBAR_OMEGA_EV * \
            np.sqrt(np.abs(self.eigenvalues))

    @property
    def wavenumbers(self) -> np.ndarray:
        """Frequencies in cm^-1 (negative: imaginary)."""
        return self.energies * EV2INVCM

    @property
    def imaginary(self) -> np.ndarray:
        """Indices of imaginary (negative curvature) modes."""
        return np.flatnonzero(self.eigenvalues < 0.0)

    def cartesian_modes(self) -> np.ndarray:
        """Normalized Cartesian displacement patterns (n_vib, N, 3)."""
        inv_sqm = np.repeat(1.0 / np.sqrt(self.masses), 3)
        c = self.modes * inv_sqm[:, None]
        c /= np.linalg.norm(c, axis=0)
        return c.T.reshape(self.nmodes, self.natoms, 3)

    def atoms(self) -> Atoms:
        """The reference geometry as ASE atoms."""
        a = Atoms(self.symbols, positions=self.positions)
        a.set_masses(self.masses)
        return a

    def summary(self) -> str:
        """Short text table of the frequencies."""
        wn = self.wavenumbers
        lines = [f"{self.nmodes} vibrational modes, "
                 f"{len(self.imaginary)} imaginary"]
        for i in range(0, self.nmodes, 8):
            lines.append(" ".join(f"{x:9.2f}" for x in wn[i:i + 8]))
        return "\n".join(lines)

    # ------------------------------------------------------------------
    def save(self, filename: Pathlike) -> None:
        """Write geometry, masses and Hessian to an npz file."""
        np.savez(filename, symbols=np.array(self.symbols),
                 positions=self.positions, masses=self.masses,
                 hessian=self.hessian)

    @classmethod
    def load(cls, filename: Pathlike) -> "NormalModes":
        """Read a file written by save() and redo the analysis."""
        with np.load(filename) as f:
            atoms = Atoms(list(f["symbols"]), positions=f["positions"])
            return cls.from_hessian(atoms, f["hessian"], masses=f["masses"])


# ----------------------------------------------------------------------
# readers for other codes (migration of existing data)
# ----------------------------------------------------------------------
_HARTREE_PER_BOHR = units.Hartree / units.Bohr
_HARTREE_PER_BOHR2 = units.Hartree / units.Bohr**2


def read_orca_hessian(filename: Pathlike) -> np.ndarray:
    """Cartesian Hessian from an ORCA .hess file, in eV/Angstrom^2."""
    lines = Path(filename).read_text().splitlines()
    start = [i for i, line in enumerate(lines)
             if line.strip() == "$hessian"][0]
    n = int(lines[start + 1].split()[0])
    h = np.zeros((n, n))
    j = start + 2
    done = 0
    while done < n:
        cols = [int(c) for c in lines[j].split()]
        j += 1
        for r in range(n):
            vals = lines[j + r].split()
            h[r, cols] = [float(v) for v in vals[1:1 + len(cols)]]
        j += n
        done = cols[-1] + 1
    return h * _HARTREE_PER_BOHR2


def read_orca_engrad(filename: Pathlike) -> np.ndarray:
    """Gradient (N, 3) from an ORCA-style .engrad file, in eV/Angstrom.

    Also reads the .engrad files written by the legacy votca python tools
    (energy "n/a").
    """
    lines = [line.strip()
             for line in Path(filename).read_text().splitlines()]
    data = [line for line in lines if line and not line.startswith("#")]
    natoms = int(data[0])
    g = np.array([float(x) for x in data[2:2 + 3 * natoms]])
    return g.reshape(natoms, 3) * _HARTREE_PER_BOHR
