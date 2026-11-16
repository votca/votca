"""Thermal vibronic line shapes and ensemble spectra.

Line shape of a transition between two displaced harmonic surfaces with
equal frequencies (Lax/Kubo generating function), the initial state in
thermal equilibrium at temperature T:

    F(E) = 1/(2 pi) Int dtau exp(i (E - E_00) tau) C(tau),     tau = t/hbar
    C(tau) = exp( sum_nu S_nu [ (n_nu + 1) e^{-i w_nu tau}
                                 + n_nu e^{+i w_nu tau} - (2 n_nu + 1) ] )
             * exp(-i lambda_cl tau - lambda_cl kT tau^2)          (classical)
             * exp(-sigma^2 tau^2 / 2) * exp(-gamma |tau|)         (broadening)

with w_nu = hbar omega_nu (eV), n_nu the Bose occupation and E_00 the
adiabatic (0-0) energy. The classical term is the high-temperature limit
of soft modes or of an outer-sphere (environment) reorganization energy
lambda_cl: a Marcus Gaussian of variance 2 lambda_cl kT. sigma (Gaussian)
and gamma (Lorentzian half width) are phenomenological broadenings, e.g.
the instrument resolution.

F is normalized, Int F dE = 1, and its mean is the vertical energy
E_vert = E_00 + lambda_tot. LineShape stores F on an energy axis relative
to E_vert, which is the natural reference for ensembles of vertical
energies (QM/MM site energies, GW quasiparticle levels in a morphology).

The integral is evaluated by one FFT (no quadrature); the time grid
follows from the energy step and the window.

The process fixes the energy axis and with it the side of E_00 on which
the vibronic progression lies. The initial state is always the one whose
modes are used (at its minimum); the sign of the gradient does not matter
(lambda_nu ~ q_nu^2).

    =================  ==================================  ============
    process            energy axis (site energies)         progression
    =================  ==================================  ============
    absorption         photon energy, S0 -> Sn (Omega_n)   above E_00
    emission           photon energy, Sn -> S0 (Omega_n)   below E_00
    electron-addition  QP level of the added electron      above E_00
                       (LUMO, ... ; inverse photoemission)
    electron-removal   QP level of the removed electron    below E_00
                       (HOMO, ... ; photoemission)
    ionization         ionization energy IP = -E_HOMO      above E_00
    electron-affinity  electron affinity EA = -E_LUMO      below E_00
    =================  ==================================  ============

The gradient is that of Omega_n for absorption/emission, of the addition
level (dE_LUMO) for electron-addition/electron-affinity and of the removal
level (dE_HOMO) for electron-removal/ionization. For emission the modes
should strictly be those of the excited state.
"""

from __future__ import annotations

import dataclasses
import warnings
from typing import Optional, Sequence

import numpy as np

from .coupling import VibronicCouplings
from .modes import KB_EV

#: processes and whether their vibronic progression lies above E_00
PROCESSES = {"absorption": True, "emission": False,
             "electron-addition": True, "electron-removal": False,
             "ionization": True, "electron-affinity": False}


def _check_process(process: str) -> bool:
    if process not in PROCESSES:
        raise ValueError(f"unknown process {process!r}; choose one of "
                         f"{', '.join(PROCESSES)}")
    return PROCESSES[process]


@dataclasses.dataclass
class LineShape:
    """Normalized line shape on an energy axis relative to E_vert.

    Attributes:
        energy: E - E_vert (eV), uniform grid
        intensity: F(E) (1/eV), Int F dE = 1
        reorganization: lambda_tot = |E_vert - E_00| (eV), quantum +
            classical modes
        temperature: K
        process: see PROCESSES and the module docstring
    """

    energy: np.ndarray
    intensity: np.ndarray
    reorganization: float
    temperature: float
    process: str = "absorption"

    @property
    def step(self) -> float:
        return float(self.energy[1] - self.energy[0])

    @property
    def zero_zero(self) -> float:
        """E_00 - E_vert (eV)."""
        above = _check_process(self.process)
        return -self.reorganization if above else self.reorganization

    def __call__(self, e: np.ndarray) -> np.ndarray:
        """F at energies relative to E_vert (linear interpolation, 0
        outside the grid)."""
        return np.interp(e, self.energy, self.intensity, left=0.0,
                         right=0.0)

    def moments(self):
        """(mean, standard deviation) relative to E_vert."""
        w = self.intensity * self.step
        mean = float((w * self.energy).sum())
        var = float((w * (self.energy - mean) ** 2).sum())
        return mean, np.sqrt(max(var, 0.0))


def _bose(w: np.ndarray, kt: float) -> np.ndarray:
    if kt <= 0.0:
        return np.zeros_like(w)
    return 1.0 / np.expm1(w / kt)


def vibronic_lineshape(couplings: Optional[VibronicCouplings],
                       temperature: float,
                       process: str,
                       lambda_classical: float = 0.0,
                       sigma: float = 0.0,
                       gamma: float = 0.0,
                       classical_below: Optional[float] = None,
                       step: float = 1e-3,
                       window: Optional[float] = None,
                       tolerance: float = 1e-6) -> LineShape:
    """Thermal Franck-Condon line shape by FFT.

    Args:
        couplings: quantum modes (None: classical/broadening only)
        temperature: K (0 allowed for the quantum part)
        process: absorption, emission, electron-addition,
            electron-removal, ionization or electron-affinity; sets the
            energy axis and the side of E_00 on which the progression lies
            (module docstring)
        lambda_classical: classical reorganization energy (eV), e.g. outer
            sphere; needs temperature > 0
        sigma: Gaussian broadening, standard deviation (eV)
        gamma: Lorentzian broadening, half width at half maximum (eV)
        classical_below: modes with hbar omega below this (eV) are treated
            classically (their lambda is added to lambda_classical)
        step: energy step (eV)
        window: total width of the energy axis (eV), centred on E_vert;
            default from the width of the spectrum
        tolerance: warn if |C(tau)| at the end of the time grid exceeds
            this (spectrum not converged: add broadening or a classical
            part, or reduce the step)

    Returns:
        LineShape relative to E_vert.
    """
    above = _check_process(process)
    kt = KB_EV * temperature
    if couplings is None:
        w = np.zeros(0)
        lam = np.zeros(0)
    else:
        if classical_below is not None:
            low, couplings = couplings.split(classical_below)
            lambda_classical += low.total_reorganization
        w = np.asarray(couplings.hbar_omega, dtype=float)
        lam = np.asarray(couplings.reorganization, dtype=float)
    if (w <= 0).any():
        raise ValueError("vibrational quanta must be positive")
    if lambda_classical > 0.0 and kt <= 0.0:
        raise ValueError("a classical reorganization energy needs T > 0")
    s = lam / w if w.size else w
    n = _bose(w, kt)
    lam_tot = float(lam.sum()) + lambda_classical

    if lam.size == 0 and lambda_classical == 0 and sigma == 0 and gamma == 0:
        raise ValueError("nothing to broaden: no modes and no broadening")

    # width: variance of the line shape (second cumulant) + broadenings
    var = float((s * w**2 * (2 * n + 1)).sum()) + 2 * lambda_classical * kt \
        + sigma**2
    width = np.sqrt(var) + 10 * gamma
    if window is None:
        window = max(2.0 * (lam_tot + 12.0 * width), 0.5)
    npts = int(2 ** np.ceil(np.log2(window / step)))
    de = window / npts

    # tau = t/hbar in 1/eV, in FFT order
    tau = 2.0 * np.pi * np.fft.fftfreq(npts, d=de)
    # log C relative to E_vert: shift by +lambda_tot (mean of F is E_vert)
    logc = 1j * lam_tot * tau
    if w.size:
        wt = np.outer(tau, w)
        logc += ((n + 1) * np.exp(-1j * wt) + n * np.exp(1j * wt)
                 - (2 * n + 1)) @ s
    logc += -1j * lambda_classical * tau - lambda_classical * kt * tau**2
    logc += -0.5 * sigma**2 * tau**2 - gamma * np.abs(tau)
    c = np.exp(logc)
    tail = np.abs(c[npts // 2])
    if tail > tolerance:
        warnings.warn(f"line shape not converged (|C| = {tail:.1e} at the "
                      "end of the time grid): add broadening, a classical "
                      "part or use a smaller step")
    f = np.fft.fftshift(np.fft.ifft(c).real) / de
    e = np.fft.fftshift(np.fft.fftfreq(npts, d=1.0 / (npts * de)))
    if not above:
        # F(E_vert - x): reverse the axis; keep the grid symmetric
        f = np.roll(f[::-1], 1)
    return LineShape(energy=e, intensity=f, reorganization=lam_tot,
                     temperature=temperature, process=process)


def ensemble_spectrum(vertical_energies: Sequence[float],
                      shape: LineShape,
                      energies: Optional[np.ndarray] = None,
                      weights: Optional[Sequence[float]] = None,
                      step: Optional[float] = None):
    """Sum of line shapes centred on vertical energies.

    I(E) = sum_i w_i F(E - E_i). With weights None every site counts 1, so
    Int I dE is the number of sites. Weights can model e.g. the surface
    sensitivity of photoemission, w_i = exp(-depth_i / IMFP).

    Args:
        vertical_energies: E_vert of each site/state (eV)
        shape: LineShape (relative to E_vert)
        energies: output grid (eV); default covers all sites
        weights: optional weight per site
        step: grid step of the default grid (default: shape.step)

    Returns:
        (energies, intensity)
    """
    ev = np.asarray(vertical_energies, dtype=float)
    wts = np.ones_like(ev) if weights is None else \
        np.asarray(weights, dtype=float)
    if energies is None:
        st = shape.step if step is None else step
        mask = shape.intensity > 1e-8 * shape.intensity.max()
        lo, hi = shape.energy[mask][[0, -1]]
        energies = np.arange(ev.min() + lo, ev.max() + hi + st, st)
    energies = np.asarray(energies, dtype=float)
    # histogram the sites on the line-shape grid, then convolve: exact on
    # the grid and O(M log M) instead of O(M N_sites)
    de = shape.step
    lo = energies.min() - shape.energy[-1] - de
    hi = energies.max() - shape.energy[0] + de
    nbins = int(np.ceil((hi - lo) / de)) + 1
    grid = lo + de * np.arange(nbins)
    pos = (ev - lo) / de
    i0 = np.floor(pos).astype(int)
    frac = pos - i0
    hist = np.zeros(nbins)
    inside = (i0 >= 0) & (i0 + 1 < nbins)
    np.add.at(hist, i0[inside], wts[inside] * (1 - frac[inside]))
    np.add.at(hist, i0[inside] + 1, wts[inside] * frac[inside])
    nconv = hist.size + shape.intensity.size - 1
    conv = np.fft.irfft(np.fft.rfft(hist, nconv) *
                        np.fft.rfft(shape.intensity, nconv), nconv)
    # conv[k] corresponds to grid[0] + shape.energy[0] + k*de
    axis = grid[0] + shape.energy[0] + de * np.arange(conv.size)
    return energies, np.interp(energies, axis, conv, left=0.0, right=0.0)
