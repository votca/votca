"""Mode-resolved vibronic couplings (displaced harmonic oscillators).

Model: the initial state has the harmonic surface of ``NormalModes`` with
its minimum at the reference geometry; the final state has the same
frequencies, displaced by the linear term g = grad(E_final - E_initial)
at the reference geometry. In mass-weighted normal coordinates Q_nu

    E_final - E_initial = E_vert + sum_nu q_nu Q_nu,
    q_nu = e_nu . (M^-1/2 g),

so the final-state minimum is displaced by dQ_nu = -q_nu / omega_nu^2 and

    lambda_nu = q_nu^2 / (2 omega_nu^2)          (reorganization energy)
    S_nu      = lambda_nu / (hbar omega_nu)       (Huang-Rhys factor).

lambda = sum_nu lambda_nu is the relaxation energy of the final state from
the vertical point (the one-sided reorganization energy; for charge
transfer lambda_ct = lambda(charging) + lambda(discharging)).

The gradient is typically
- an excitation energy (BSE Omega_n, from pyxtp ``xtp(excited_state=...)``
  forces minus ground-state forces),
- a quasiparticle energy (dE_QP/dR: electron addition to the LUMO, removal
  from the HOMO),
- or any total-energy difference.
"""

from __future__ import annotations

import dataclasses
import warnings
from typing import Optional

import numpy as np
from ase import Atoms

from .modes import EV2INVCM, NormalModes


@dataclasses.dataclass
class VibronicCouplings:
    """Huang-Rhys factors and reorganization energies per mode.

    Attributes
    ----------
    hbar_omega
        vibrational quanta (n), eV; only real modes are kept
    reorganization
        lambda_nu (n), eV
    mode_index
        index of each entry in the NormalModes it came from
    dropped_lambda
        reorganization energy of modes that were left out
        (imaginary modes, modes below a cutoff), eV
    net_force
        norm of the net force/torque part of the gradient that
        is not carried by vibrations (diagnostic), eV/Angstrom

    """

    hbar_omega: np.ndarray
    reorganization: np.ndarray
    mode_index: np.ndarray
    dropped_lambda: float = 0.0
    net_force: float = 0.0

    @property
    def huang_rhys(self) -> np.ndarray:
        """S_nu = lambda_nu / hbar omega_nu."""
        return self.reorganization / self.hbar_omega

    @property
    def total_reorganization(self) -> float:
        """Total reorganization energy, sum_nu lambda_nu (eV)."""
        return float(self.reorganization.sum())

    @property
    def wavenumbers(self) -> np.ndarray:
        """Mode frequencies in cm^-1."""
        return self.hbar_omega * EV2INVCM

    def split(self, hbar_omega_max: float):
        """(low, high): couplings of modes below / at or above a quantum.

        Used to treat soft modes classically (see lineshape).
        """
        low = self.hbar_omega < hbar_omega_max

        def take(mask):
            return VibronicCouplings(self.hbar_omega[mask],
                                     self.reorganization[mask],
                                     self.mode_index[mask])
        return take(low), take(~low)

    def summary(self, top: int = 10) -> str:
        """Text table of the strongest modes."""
        order = np.argsort(self.reorganization)[::-1][:top]
        lines = [f"lambda = {self.total_reorganization:.4f} eV, "
                 f"sum S = {self.huang_rhys.sum():.3f}, "
                 f"{len(self.hbar_omega)} modes"
                 + (f", {self.dropped_lambda:.4f} eV dropped"
                    if self.dropped_lambda else ""),
                 "   cm^-1   lambda/meV       S"]
        for k in order:
            lines.append(f"{self.wavenumbers[k]:8.1f} "
                         f"{1000 * self.reorganization[k]:12.2f} "
                         f"{self.huang_rhys[k]:7.3f}")
        return "\n".join(lines)


def couplings_from_gradient(modes: NormalModes, gradient: np.ndarray,
                            min_wavenumber: float = 0.0,
                            ) -> VibronicCouplings:
    """Couplings from the gradient of the energy gap.

    Parameters
    ----------
    modes
        normal modes of the initial state at its minimum
    gradient
        grad(E_final - E_initial) at modes.positions, (N, 3) or
        (3N,), eV/Angstrom. With forces F: gradient = -(F_final -
        F_initial); see couplings_from_forces.
    min_wavenumber
        modes below this frequency (cm^-1) are left out
        (their lambda is reported as dropped_lambda). Imaginary modes
        are always left out with a warning.

    """
    g = np.asarray(gradient, dtype=float).ravel()
    if g.size != 3 * modes.natoms:
        raise ValueError(f"gradient has {g.size} components, expected "
                         f"{3 * modes.natoms}")
    gm = g / np.repeat(np.sqrt(modes.masses), 3)
    q = modes.modes.T @ gm
    rest = gm - modes.modes @ q
    net = float(np.linalg.norm(rest * np.repeat(np.sqrt(modes.masses), 3)))

    k = modes.eigenvalues
    real = k > 0.0
    lam = np.zeros_like(q)
    lam[real] = q[real] ** 2 / (2.0 * k[real])
    if (~real).any():
        warnings.warn(f"{(~real).sum()} imaginary mode(s) left out of the "
                      "couplings; the reference is not a minimum")
    keep = real & (modes.wavenumbers >= min_wavenumber)
    dropped = float(lam[real & ~keep].sum())
    idx = np.flatnonzero(keep)
    return VibronicCouplings(hbar_omega=modes.energies[idx],
                             reorganization=lam[idx], mode_index=idx,
                             dropped_lambda=dropped, net_force=net)


def couplings_from_forces(modes: NormalModes, final_forces: np.ndarray,
                          initial_forces: Optional[np.ndarray] = None,
                          **kwargs) -> VibronicCouplings:
    """Couplings from the forces (eV/Angstrom) of both states.

    The initial-state forces should vanish at a minimum; passing them
    removes the effect of a slightly unconverged geometry on the gap
    gradient.
    """
    g = -np.asarray(final_forces, dtype=float)
    if initial_forces is not None:
        g = g + np.asarray(initial_forces, dtype=float)
    return couplings_from_gradient(modes, g, **kwargs)


def read_gap_gradient(filename) -> np.ndarray:
    """Gradient of an energy gap (N, 3) in eV/Angstrom from a file.

    - ``.orb`` (dftgwbse with ``excited_state_forces``): grad Omega =
      F_ground - F_excited, from the datasets ``forces`` and
      ``excited_state_forces`` (Hartree/Bohr)
    - ``.engrad``: ORCA-style gradient (Hartree/Bohr), e.g. QP-level
      gradients of the legacy votca python tools
    - ``.npy`` or text: (N, 3) array in eV/Angstrom
    """
    from pathlib import Path

    from ase import units

    from .modes import read_orca_engrad

    path = Path(filename)
    if path.suffix == ".orb":
        import h5py
        with h5py.File(path, "r") as f:
            orb = f["QMdata"] if "QMdata" in f else f
            if "excited_state_forces" not in orb:
                raise ValueError(f"{path} has no excited_state_forces; run "
                                 "dftgwbse with excited_state_forces set")
            g = orb["forces"][()] - orb["excited_state_forces"][()]
        return g * units.Hartree / units.Bohr
    if path.suffix == ".engrad":
        return read_orca_engrad(path)
    if path.suffix == ".npy":
        return np.load(path).reshape(-1, 3)
    return np.loadtxt(path).reshape(-1, 3)


def _kabsch(ref: np.ndarray, mov: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Rotate and translate mov onto ref (mass-weighted best fit)."""
    w = w / w.sum()
    rc = (w[:, None] * ref).sum(0)
    mc = (w[:, None] * mov).sum(0)
    a = (w[:, None] * (mov - mc)).T @ (ref - rc)
    u, _, vt = np.linalg.svd(a)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    return (mov - mc) @ r + rc


def couplings_from_geometry(modes: NormalModes, final: Atoms,
                            min_wavenumber: float = 0.0,
                            ) -> VibronicCouplings:
    """Couplings from the optimized geometry of the final state.

    dQ = modes^T M^1/2 (R_final - R_initial) after a mass-weighted fit;
    lambda_nu = omega_nu^2 dQ_nu^2 / 2. Agrees with the gradient route if
    the final state is harmonic with the initial-state frequencies; the
    difference measures anharmonicity and frequency changes.
    """
    r = _kabsch(modes.positions, final.get_positions(), modes.masses)
    dx = (r - modes.positions).ravel() * np.repeat(np.sqrt(modes.masses), 3)
    dq = modes.modes.T @ dx
    k = modes.eigenvalues
    real = k > 0.0
    lam = np.where(real, 0.5 * k * dq**2, 0.0)
    keep = real & (modes.wavenumbers >= min_wavenumber)
    idx = np.flatnonzero(keep)
    return VibronicCouplings(hbar_omega=modes.energies[idx],
                             reorganization=lam[idx], mode_index=idx,
                             dropped_lambda=float(lam[real & ~keep].sum()))
