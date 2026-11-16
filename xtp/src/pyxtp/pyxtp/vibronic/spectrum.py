"""Vibronically resolved spectra of several levels (UPS, IPES).

``QPLevels.from_orb`` reads the quasiparticle energies and their analytic
gradients (dftgwbse option ``qp_gradients``, e.g. ``homo-5:lumo+5``) of a
molecule. ``level_spectrum`` turns them into a spectrum without manual
steps: for every level that fits the process (occupied levels for
electron-removal/ionization, empty levels for electron-addition/
electron-affinity) the couplings and the thermal line shape are computed
and placed at the vertical energies of the level -- the QP energy of the
molecule, or the site energies of a morphology -- and the contributions
are added.

All levels count equally (one state per level and site): photoemission
cross sections and their dependence on the photon energy are not included.
"""

from __future__ import annotations

import dataclasses
import re
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Union

import numpy as np
from ase import units

from .coupling import VibronicCouplings, couplings_from_gradient
from .lineshape import (PROCESSES, _check_process, ensemble_spectrum,
                        vibronic_lineshape)
from .modes import NormalModes

Pathlike = Union[str, Path]

#: processes of occupied (True) and empty (False) levels, and whether the
#: energy axis is minus the level energy (IP, EA)
_OCCUPIED = {"electron-removal": True, "ionization": True,
             "electron-addition": False, "electron-affinity": False}
_NEGATED = {"ionization", "electron-affinity"}


def level_label(level: int, homo: int) -> str:
    """HOMO-2, HOMO, LUMO, LUMO+3, ..."""
    if level <= homo:
        return "HOMO" if level == homo else f"HOMO-{homo - level}"
    return "LUMO" if level == homo + 1 else f"LUMO+{level - homo - 1}"


def parse_level(text: str, homo: int) -> int:
    """Absolute MO index of 'homo-2', 'lumo+3', 'HOMO' or '17'."""
    t = text.strip().lower().replace(" ", "")
    for name, base in (("homo", homo), ("lumo", homo + 1)):
        if t.startswith(name):
            rest = t[len(name):]
            return base + (int(rest) if rest else 0)
    return int(t)


@dataclasses.dataclass
class QPLevels:
    """QP levels of a molecule with the gradients of their energies.

    Attributes:
        levels: absolute MO indices (0-based)
        homo: index of the HOMO
        energies: QP energies (eV), one per level
        gradients: dE/dR (n_levels, N, 3), eV/Angstrom
    """

    levels: np.ndarray
    homo: int
    energies: np.ndarray
    gradients: np.ndarray

    @classmethod
    def from_orb(cls, filename: Pathlike) -> "QPLevels":
        """Read an orb file of dftgwbse with ``qp_gradients``."""
        import h5py
        with h5py.File(filename, "r") as f:
            orb = f["QMdata"]
            if "qp_gradients" not in orb:
                raise ValueError(f"{filename} has no QP gradients; run "
                                 "dftgwbse with qp_gradients, e.g. "
                                 "homo-3:lumo+3")
            levels = np.array(orb["qp_gradient_levels"][()],
                              dtype=int).ravel()
            g = np.array(orb["qp_gradients"][()], dtype=float)
            homo = int(orb.attrs["occupied_levels"][0]) - 1
            qpmin = int(orb.attrs["qpmin"][0])
            qp = np.array(orb["QPpert_energies"][()], dtype=float).ravel()
        g = g.reshape(len(levels), -1, 3) * units.Hartree / units.Bohr
        return cls(levels=levels, homo=homo,
                   energies=qp[levels - qpmin] * units.Hartree, gradients=g)

    @property
    def labels(self) -> List[str]:
        return [level_label(int(lv), self.homo) for lv in self.levels]

    def select(self, occupied: bool) -> "QPLevels":
        """The occupied or the empty levels."""
        mask = (self.levels <= self.homo) == occupied
        return QPLevels(self.levels[mask], self.homo, self.energies[mask],
                        self.gradients[mask])

    def couplings(self, modes: NormalModes, min_wavenumber: float = 0.0
                  ) -> List[VibronicCouplings]:
        """Vibronic couplings of every level."""
        return [couplings_from_gradient(modes, g,
                                        min_wavenumber=min_wavenumber)
                for g in self.gradients]


@dataclasses.dataclass
class Spectrum:
    """A spectrum with the contributions of the individual levels.

    Attributes:
        energy: energy axis (eV)
        contributions: label -> intensity on the energy axis
        process: physical process (sets the meaning of the axis)
        reorganization: label -> lambda (eV)
        notes: remarks on how each contribution was placed
    """

    energy: np.ndarray
    contributions: Dict[str, np.ndarray]
    process: str
    reorganization: Dict[str, float] = dataclasses.field(default_factory=dict)
    notes: Dict[str, str] = dataclasses.field(default_factory=dict)

    @property
    def total(self) -> np.ndarray:
        return np.sum(list(self.contributions.values()), axis=0)

    @property
    def axis_label(self) -> str:
        return {"absorption": "Photon energy (eV)",
                "emission": "Photon energy (eV)",
                "electron-addition": "Energy (eV)",
                "electron-removal": "Energy (eV)",
                "ionization": "Ionization energy (eV)",
                "electron-affinity": "Electron affinity (eV)"}[self.process]

    def broadened(self, sigma: float) -> "Spectrum":
        """Convolution with a Gaussian of standard deviation sigma (eV),
        e.g. the instrument resolution."""
        if sigma <= 0.0:
            return self
        de = self.energy[1] - self.energy[0]
        half = int(np.ceil(6.0 * sigma / de))
        x = de * np.arange(-half, half + 1)
        kernel = np.exp(-0.5 * (x / sigma) ** 2)
        kernel /= kernel.sum()
        # pad so that intensity near the ends of the axis is not folded
        pad = half
        e = np.concatenate([self.energy[0] + de * np.arange(-pad, 0),
                            self.energy,
                            self.energy[-1] + de * np.arange(1, pad + 1)])
        out = {k: np.convolve(np.pad(v, pad), kernel, mode="same")
               for k, v in self.contributions.items()}
        return Spectrum(e, out, self.process, dict(self.reorganization),
                        dict(self.notes))

    def onset(self, side: Optional[str] = None,
              intensity: Optional[np.ndarray] = None):
        """Onset by the tangent construction used for UPS/IPES.

        The tangent at the steepest point of the edge of the outermost
        peak is extrapolated to zero intensity.

        Args:
            side: "low" (edge on the low-energy side of the first peak;
                default for electron-addition, ionization, absorption) or
                "high" (high-energy edge of the last peak; default for
                electron-removal, electron-affinity, emission)
            intensity: curve to analyse (default: the total)

        Returns:
            dict with onset, peak, the tangent point (x0, y0) and slope
        """
        y = self.total if intensity is None else np.asarray(intensity)
        e = self.energy
        if side is None:
            side = "low" if self.process in ("electron-addition",
                                             "ionization",
                                             "absorption") else "high"
        if side not in ("low", "high"):
            raise ValueError("side must be 'low' or 'high'")
        if side == "high":
            res = Spectrum(-e[::-1], {"t": y[::-1]}, "ionization").onset(
                "low")
            return {"onset": -res["onset"], "peak": -res["peak"],
                    "x0": -res["x0"], "y0": res["y0"],
                    "slope": -res["slope"]}
        # first local maximum that is a real peak (> 5 % of the maximum)
        ymax = y.max()
        peak = None
        for i in range(1, len(y) - 1):
            if y[i] >= y[i - 1] and y[i] > y[i + 1] and y[i] > 0.05 * ymax:
                peak = i
                break
        if peak is None:
            peak = int(np.argmax(y))
        dy = np.gradient(y, e)
        i0 = int(np.argmax(dy[:peak + 1]))
        slope = float(dy[i0])
        return {"onset": float(e[i0] - y[i0] / slope), "peak": float(e[peak]),
                "x0": float(e[i0]), "y0": float(y[i0]), "slope": slope}

    def save(self, filename: Pathlike, header: str = "") -> None:
        """Text file: energy, total, one column per contribution."""
        labels = list(self.contributions)
        cols = [self.energy, self.total] + [self.contributions[k]
                                            for k in labels]
        lines = [header] if header else []
        lines.append(f"process {self.process}")
        for k in labels:
            note = self.notes.get(k, "")
            lam = self.reorganization.get(k)
            lines.append(f"{k}: " + (f"lambda {lam:.4f} eV" if lam is not None
                                     else "") + (f"; {note}" if note else ""))
        lines.append("E/eV total " + " ".join(labels))
        np.savetxt(filename, np.column_stack(cols), header="\n".join(lines))

    @classmethod
    def load(cls, filename: Pathlike) -> "Spectrum":
        """Read a file written by save()."""
        header = []
        with open(filename) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                header.append(line[1:].strip())
        data = np.loadtxt(filename, ndmin=2)
        found = re.search(r"process (\S+)", " ".join(header))
        process = found.group(1) if found else "absorption"
        if data.shape[1] == 2:  # output of 'xtp_vibronic lineshape'
            return cls(data[:, 0], {"Spectrum": data[:, 1]}, process)
        labels = header[-1].split()[2:]
        if len(labels) != data.shape[1] - 2:
            labels = [f"c{i}" for i in range(data.shape[1] - 2)]
        lam = {}
        for h in header:
            for k in labels:
                if h.startswith(k + ": lambda "):
                    lam[k] = float(h.split()[2])
        return cls(data[:, 0], {k: data[:, 2 + i]
                                for i, k in enumerate(labels)}, process, lam)


def level_spectrum(modes: NormalModes, levels: QPLevels, process: str,
                   temperature: float,
                   site_energies: Optional[Mapping[Union[int, str],
                                                   Sequence[float]]] = None,
                   weights: Optional[Sequence[float]] = None,
                   select: Optional[Sequence[Union[int, str]]] = None,
                   lambda_classical: float = 0.0, sigma: float = 0.0,
                   gamma: float = 0.0,
                   classical_below: Optional[float] = None,
                   min_wavenumber: float = 0.0,
                   energies: Optional[np.ndarray] = None,
                   step: float = 2e-3) -> Spectrum:
    """Vibronic spectrum of the QP levels of a molecule or a morphology.

    Args:
        modes: normal modes of the neutral ground state
        levels: QP levels with gradients (QPLevels.from_orb)
        process: electron-removal or ionization (occupied levels; axis
            E_QP or -E_QP), electron-addition or electron-affinity (empty
            levels)
        temperature: K
        site_energies: vertical energies of the sites of a morphology per
            level (key: MO index or a label such as 'homo-1', 'LUMO'), on
            the energy axis of the process. Levels without an entry use
            the distribution of the nearest level that has one, shifted
            rigidly by the difference of the molecular QP energies.
            None: one molecule, at its QP energies.
        weights: weight per site (e.g. exp(-depth/IMFP)), the same for all
            levels
        select: use only these levels (indices or labels)
        lambda_classical, sigma, gamma, classical_below: see
            vibronic_lineshape (intrinsic line shape of every level)
        min_wavenumber: leave out modes below this (cm^-1)
        energies: output energy axis (default covers all contributions)
        step: energy step of the line shapes and the default axis (eV)

    Returns:
        Spectrum with one contribution per level.
    """
    _check_process(process)
    if process not in _OCCUPIED:
        raise ValueError(f"process {process} has no QP levels; use "
                         f"{', '.join(_OCCUPIED)}")
    sub = levels.select(_OCCUPIED[process])
    if select is not None:
        wanted = {parse_level(str(s), levels.homo) for s in select}
        mask = np.isin(sub.levels, list(wanted))
        sub = QPLevels(sub.levels[mask], sub.homo, sub.energies[mask],
                       sub.gradients[mask])
    if len(sub.levels) == 0:
        raise ValueError(f"no {'occupied' if _OCCUPIED[process] else 'empty'}"
                         f" levels with gradients for {process}")
    sign = -1.0 if process in _NEGATED else 1.0
    vertical = sign * sub.energies

    sites: Dict[int, np.ndarray] = {}
    for key, val in (site_energies or {}).items():
        sites[parse_level(str(key), levels.homo)] = np.asarray(val,
                                                               dtype=float)
    unknown = set(sites) - set(int(lv) for lv in sub.levels)
    if unknown:
        raise ValueError("site energies given for levels without gradients "
                         f"or of the wrong kind for {process}: "
                         f"{sorted(unknown)}")

    shapes = []
    lam = {}
    for g in sub.gradients:
        c = couplings_from_gradient(modes, g, min_wavenumber=min_wavenumber)
        shapes.append(vibronic_lineshape(
            c, temperature, process, lambda_classical=lambda_classical,
            sigma=sigma, gamma=gamma, classical_below=classical_below,
            step=step))
    placed = []
    notes = {}
    for i, lv in enumerate(int(x) for x in sub.levels):
        label = level_label(lv, sub.homo)
        lam[label] = shapes[i].reorganization
        if not sites:
            placed.append(np.array([vertical[i]]))
            notes[label] = f"molecular QP energy {sub.energies[i]:.4f} eV"
        elif lv in sites:
            placed.append(sites[lv])
            notes[label] = f"{len(sites[lv])} site energies"
        else:
            ref = min(sites, key=lambda r: abs(r - lv))
            j = int(np.flatnonzero(sub.levels == ref)[0])
            shift = vertical[i] - vertical[j]
            placed.append(sites[ref] + shift)
            notes[label] = (f"site energies of {level_label(ref, sub.homo)} "
                            f"shifted by {shift:+.4f} eV (molecular QP "
                            "energy difference)")
    if energies is None:
        lo = min(p.min() + s.energy[s.intensity > 1e-8 * s.intensity.max()][0]
                 for p, s in zip(placed, shapes))
        hi = max(p.max() + s.energy[s.intensity > 1e-8 * s.intensity.max()][-1]
                 for p, s in zip(placed, shapes))
        energies = np.arange(lo, hi + step, step)
    contributions = {}
    for i, lv in enumerate(int(x) for x in sub.levels):
        w = None
        if weights is not None and len(placed[i]) == len(weights):
            w = weights
        _, inten = ensemble_spectrum(placed[i], shapes[i], energies=energies,
                                     weights=w)
        contributions[level_label(lv, sub.homo)] = inten
    return Spectrum(np.asarray(energies), contributions, process, lam, notes)


def add_gradients(orbfile: Pathlike, qp_levels: Optional[str] = None,
                  excited_state: Optional[str] = None, options=None,
                  settings: Optional[Mapping[str, str]] = None,
                  nthreads: int = 1) -> None:
    """Add gradients to the orb file of an earlier dftgwbse run.

    Runs dftgwbse with ``tasks=gradients``: no DFT and no GW-BSE, only the
    gradients, from the data in the orb file. QP levels whose gradients
    are already stored are not recomputed, so levels can be added later;
    the forces of another excited state replace the stored ones.

    Args:
        orbfile: orb file of the earlier run; it is updated in place
        qp_levels: e.g. 'homo-9:lumo+9' (dftgwbse option qp_gradients)
        excited_state: e.g. 's2' (dftgwbse option excited_state_forces;
            needs the ground-state forces in the file)
        options: pyxtp Options of the earlier run (e.g. calc.options)
        settings: or its settings as for xtp.set(); the DFT and GW
            options must be those of the earlier run (basis sets and
            functional are checked, the GW setup through the QP energies
            rebuilt for the gradient)
        nthreads: threads
    """
    import os

    from pyxtp_binds import xtp_binds

    from ..capture_standard_output import capture_standard_output
    from ..xtp import xtp

    if qp_levels is None and excited_state is None:
        raise ValueError("give qp_levels and/or excited_state")
    orb = Path(orbfile).absolute()
    if not orb.exists():
        raise FileNotFoundError(orb)
    calc = xtp(nthreads=nthreads)
    if options is not None:
        calc.set(**options.to_flat_dict())
    calc.set(**dict(settings or {}))
    calc.set(tasks="gradients")
    calc.options.job_name = orb.stem
    if qp_levels is not None:
        calc.options.qp_gradients = qp_levels
    if excited_state is not None:
        calc.options.excited_state_forces = excited_state
        calc.options.dftpackage.xtpdft.compute_forces = "true"
    cwd = os.getcwd()
    os.chdir(orb.parent)
    try:
        xml = f"{orb.stem}_gradients.xml"
        calc.options._write_xml(xml)
        output = capture_standard_output(
            xtp_binds.call_tool, "dftgwbse", nthreads,
            Path(xml).absolute().as_posix())
        with open(f"{orb.stem}_gradients.out", "w") as handle:
            handle.write(output)
    finally:
        os.chdir(cwd)


__all__ = ["QPLevels", "Spectrum", "level_spectrum", "level_label",
           "parse_level", "add_gradients", "PROCESSES"]
