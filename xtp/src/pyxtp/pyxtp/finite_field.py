"""
Excited-state dipole moments and polarizabilities by finite fields.

The excitation energies Omega_n of GW-BSE (or the DFT ground-state energy)
are computed in small uniform electric fields with ``dftgwbse``
(``dftpackage.externalfield``) and differentiated numerically. This gives
the *relaxed* change of the dipole moment on excitation, including the
response of the orbitals, of the quasiparticle energies and of the
screened interaction W -- unlike the excited-state density of
``Orbitals::DensityMatrixFull``, which is unrelaxed.

Conventions: VOTCA adds -F.r to the one-electron Hamiltonian for a field
F given in ``externalfield`` (Hartree/bohr). With E(F) = E(0) + mu.F -
1/2 F.alpha.F one gets

    mu_el   = dE_DFT/dF       (electronic dipole w.r.t. the centre of mass)
    dmu_n   = dOmega_n/dF     (change of the dipole on excitation)
    dalpha_n,ii = -d2 Omega_n/dF_i^2

Everything is in atomic units (e bohr, bohr^3) unless stated otherwise.

States are matched between the field calculations and the zero-field
calculation by the overlap of their BSE amplitudes, after transforming
the amplitudes into the zero-field molecular orbitals. States that mix in
the field (near degeneracies) are reported as groups; the derivative of
the sum of their energies is reliable even when the individual ones are
not.
"""

from __future__ import annotations

import contextlib
import os
# only used to write XML, never to parse it
import xml.etree.ElementTree as ET  # nosec B405
import dataclasses
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import h5py
import numpy as np

AU2DEBYE = 2.541746473
_AXES = {"x": 0, "y": 1, "z": 2}


# --------------------------------------------------------------------------
# reading the orb files
# --------------------------------------------------------------------------
@dataclasses.dataclass
class OrbData:
    """The parts of an orb (HDF5) file needed here (atomic units)."""

    dft_energy: float
    omega: np.ndarray  # (nstates,)
    X: np.ndarray  # (nstates, nv, nc) resonant amplitudes
    Y: Optional[np.ndarray]  # anti-resonant amplitudes, None for TDA
    mos: np.ndarray  # (nao, nmo), columns are MOs
    vmin: int
    homo: int
    cmax: int
    elements: List[str]
    positions: np.ndarray  # (natoms, 3), bohr
    charges: np.ndarray  # nuclear charges minus ECP core charges


def read_orb(filename: os.PathLike, spin: str = "singlet") -> OrbData:
    """Read energies, BSE amplitudes and MOs from an orb file."""
    with h5py.File(filename, "r") as h5:
        q = h5["QMdata"]
        homo = int(q.attrs["occupied_levels"][0]) - 1
        vmin = int(q.attrs["bse_vmin"][0])
        cmax = int(q.attrs["bse_cmax"][0])
        tda = bool(q.attrs["useTDA"][0])
        nv, nc = homo - vmin + 1, cmax - homo
        group = q[f"BSE_{spin}"]
        omega = group["eigenvalues"][()].flatten()
        nstates = omega.size

        def amplitudes(name):
            raw = group[name][()]
            if raw.size == 0 or raw.shape[0] != nv * nc:
                return None
            # VOTCA's compound index of a (v, c) pair is v * nc + c
            return raw.T.reshape(nstates, nv, nc)

        X = amplitudes("eigenvectors")
        if X is None:
            raise ValueError(f"{filename}: no BSE {spin} amplitudes")
        Y = None if tda else amplitudes("eigenvectors2")
        mos = q["mos"]["eigenvectors"][()]
        atoms = q["qmmolecule"]["qmatoms"][()].flatten()
        elements = [a["element"].decode() if isinstance(a["element"], bytes)
                    else str(a["element"]) for a in atoms]
        positions = np.array([[a["posX"], a["posY"], a["posZ"]]
                              for a in atoms], dtype=float)
        charges = np.array([a["nuccharge"] - a["ecpcharge"] for a in atoms],
                           dtype=float)
        return OrbData(float(q.attrs["qm_energy"][0]), omega, X, Y, mos,
                       vmin, homo, cmax, elements, positions, charges)


# --------------------------------------------------------------------------
# state matching
# --------------------------------------------------------------------------
def mo_rotation(c_ref: np.ndarray, c_new: np.ndarray) -> np.ndarray:
    """
    Return U with C_new = C_ref U.

    For a complete MO set the AO overlap is S = (C C^T)^-1, so that
    U = C_ref^T S C_new = C_ref^-1 C_new needs no integrals.
    """
    return np.linalg.solve(c_ref, c_new)


def state_overlaps(ref: OrbData, new: OrbData) -> np.ndarray:
    """
    Return |<n_ref|m_new>| for all pairs of states (normalized amplitudes).

    The amplitudes of `new` are expressed in the MOs of `ref` (occupied
    and virtual BSE windows) before the overlap is taken.
    """
    if (ref.vmin, ref.homo, ref.cmax) != (new.vmin, new.homo, new.cmax):
        raise ValueError("the BSE windows of the calculations differ")
    U = mo_rotation(ref.mos, new.mos)
    occ = slice(ref.vmin, ref.homo + 1)
    vir = slice(ref.homo + 1, ref.cmax + 1)
    Uoo, Uvv = U[occ, occ], U[vir, vir]

    def vectors(d: OrbData, rotate: bool) -> np.ndarray:
        blocks = [d.X] if d.Y is None else [d.X, d.Y]
        if rotate:
            blocks = [np.einsum("iv,nvc,ac->nia", Uoo, b, Uvv) for b in blocks]
        v = np.concatenate([b.reshape(b.shape[0], -1) for b in blocks], axis=1)
        return v / np.linalg.norm(v, axis=1, keepdims=True)

    return np.abs(vectors(ref, False) @ vectors(new, True).T)


def match_states(overlap: np.ndarray, nref: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Assign a state of the other calculation to each reference state.

    The first `nref` reference states are assigned one-to-one with maximum
    total overlap.

    Returns
    -------
    cols : np.ndarray
        Index of the assigned state for each reference state.
    overlaps : np.ndarray
        The overlaps of the assigned pairs.

    """
    sub = overlap[:nref]
    try:
        from scipy.optimize import linear_sum_assignment

        rows, cols = linear_sum_assignment(-sub)
        order = np.argsort(rows)
        cols = cols[order]
    except ImportError:  # greedy fallback
        cols = -np.ones(nref, dtype=int)
        taken = set()
        for flat in np.argsort(-sub, axis=None):
            n, m = divmod(int(flat), sub.shape[1])
            if cols[n] < 0 and m not in taken:
                cols[n] = m
                taken.add(m)
    return cols, sub[np.arange(nref), cols]


def mixing_groups(overlaps: Sequence[np.ndarray], nref: int,
                  threshold: float) -> List[List[int]]:
    """
    Group reference states that mix in any of the field calculations.

    Two reference states belong to the same group if a state of a field
    calculation has an overlap above (1 - threshold) with both, i.e. if
    neither of them is mapped cleanly. Isolated states form groups of one.
    """
    parent = list(range(nref))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    weak = 1.0 - threshold
    for ov in overlaps:
        sub = ov[:nref]
        for m in range(sub.shape[1]):
            partners = [n for n in range(nref) if sub[n, m] > weak]
            for n in partners[1:]:
                parent[find(n)] = find(partners[0])
    groups: Dict[int, List[int]] = {}
    for n in range(nref):
        groups.setdefault(find(n), []).append(n)
    return sorted(groups.values())


# --------------------------------------------------------------------------
# finite differences
# --------------------------------------------------------------------------
def central_first(values: Dict[int, float], h: float, order: int) -> float:
    """d/dF from values at multiples k*h of the step (k = -2..2)."""
    if order == 2:
        return (values[1] - values[-1]) / (2 * h)
    return (8 * (values[1] - values[-1]) - (values[2] - values[-2])) / (12 * h)


def central_second(values: Dict[int, float], h: float, order: int) -> float:
    """d2/dF2 from values at multiples k*h of the step (k = -2..2)."""
    if order == 2:
        return (values[1] + values[-1] - 2 * values[0]) / h ** 2
    return (-(values[2] + values[-2]) + 16 * (values[1] + values[-1])
            - 30 * values[0]) / (12 * h ** 2)


def center_of_mass(elements: Sequence[str], positions: np.ndarray) -> np.ndarray:
    """
    Return the centre of mass.

    As VOTCA's QMMolecule::getPos, this is the origin of the dipole
    operator coupled to externalfield.
    """
    from ase.data import atomic_masses, atomic_numbers

    m = np.array([atomic_masses[atomic_numbers[e]] for e in elements])
    return (m[:, None] * positions).sum(axis=0) / m.sum()


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
@dataclasses.dataclass
class FiniteFieldResult:
    """
    Results of a finite-field calculation in atomic units.

    Vectors are full 3-vectors with NaN for components that were not
    computed.
    """

    components: str
    field: float
    order: int
    spin: str
    omega: np.ndarray  # zero-field excitation energies (Hartree)
    ground_dipole: np.ndarray  # total ground-state dipole (w.r.t. COM)
    ground_alpha: np.ndarray  # diagonal of the ground-state polarizability
    delta_dipole: np.ndarray  # (nstates, 3)
    delta_alpha: np.ndarray  # (nstates, 3) diagonal
    min_overlap: np.ndarray  # (nstates,) smallest matching overlap
    groups: List[List[int]]
    group_delta_dipole: List[np.ndarray] = dataclasses.field(
        default_factory=list)
    linearity: Optional[np.ndarray] = None  # (nstates, 3), order 4 only
    overlap_threshold: float = 0.9

    @property
    def excited_dipole(self) -> np.ndarray:
        """Total excited-state dipoles, (nstates, 3)."""
        return self.ground_dipole[None, :] + self.delta_dipole

    def reliable(self, threshold: Optional[float] = None) -> np.ndarray:
        """
        Flag the states that can be trusted individually.

        These are states that are isolated and matched with an overlap of
        at least `threshold` (default: the threshold of the calculation).
        """
        if threshold is None:
            threshold = self.overlap_threshold
        single = {g[0] for g in self.groups if len(g) == 1}
        return np.array([(n in single) and (self.min_overlap[n] >= threshold)
                         for n in range(self.omega.size)])

    def to_xml(self) -> ET.Element:
        """
        Return the results as an XML tree.

        Values are in atomic units; only the computed field components
        appear.
        """
        comp = [_AXES[c] for c in self.components]

        def vec(tag, v, names=("x", "y", "z")):
            el = ET.Element(tag)
            for i in comp:
                el.set(names[i], f"{v[i]:.8f}")
            return el

        tensor = ("xx", "yy", "zz")
        root = ET.Element("finite_field", {
            "units": "atomic", "spin": self.spin,
            "field": f"{self.field:g}", "order": str(self.order),
            "components": self.components,
            "convention": "dmu = dOmega/dF with F = dftpackage.externalfield"})
        ground = ET.SubElement(root, "ground_state")
        ground.append(vec("dipole", self.ground_dipole))
        ground.append(vec("polarizability", self.ground_alpha, tensor))
        ok = self.reliable()
        for n in range(self.omega.size):
            st = ET.SubElement(root, "state", {
                "index": str(n + 1),
                "omega": f"{self.omega[n]:.8f}",
                "overlap": f"{self.min_overlap[n]:.4f}",
                "reliable": "true" if ok[n] else "false"})
            st.append(vec("delta_dipole", self.delta_dipole[n]))
            st.append(vec("dipole", self.excited_dipole[n]))
            st.append(vec("delta_polarizability", self.delta_alpha[n], tensor))
            if self.linearity is not None:
                st.append(vec("linearity", self.linearity[n]))
        for g, dg in zip(self.groups, self.group_delta_dipole):
            if len(g) > 1:
                grp = ET.SubElement(root, "mixing_group", {
                    "states": " ".join(str(n + 1) for n in g)})
                grp.append(vec("delta_dipole_sum", dg))
        return root

    def table(self) -> str:
        """Human-readable summary."""
        comp = [_AXES[c] for c in self.components]
        lines = [
            f"Finite-field excited-state dipoles ({self.spin}s, field step "
            f"{self.field:g} a.u., {'2nd' if self.order == 2 else '4th'}-order"
            " central differences)",
            "Ground-state dipole (a.u., w.r.t. centre of mass): "
            + "  ".join(f"{'xyz'[i]} {self.ground_dipole[i]:+.4f}" for i in comp)
            + f"   |mu| = {np.linalg.norm(self.ground_dipole[comp]):.4f}"
            + f" ({np.linalg.norm(self.ground_dipole[comp]) * AU2DEBYE:.3f} D)",
            "",
            "state  Omega(eV)  " + "  ".join(f"dmu_{'xyz'[i]:<5s}" for i in comp)
            + "  |dmu|(a.u.)  |dmu|(D)  " + "  ".join(
                f"dalpha_{'xyz'[i]}{'xyz'[i]}" for i in comp) + "  overlap",
        ]
        ok = self.reliable()
        for n in range(self.omega.size):
            d = self.delta_dipole[n, comp]
            flag = "" if ok[n] else "  *"
            lines.append(
                f"{n + 1:5d}  {self.omega[n] * 27.211386245988:9.4f}  "
                + "  ".join(f"{x:+9.4f}" for x in d)
                + f"  {np.linalg.norm(d):11.4f}  {np.linalg.norm(d) * AU2DEBYE:8.3f}  "
                + "  ".join(f"{self.delta_alpha[n, i]:+11.2f}" for i in comp)
                + f"  {self.min_overlap[n]:7.3f}{flag}")
        mixed = [g for g in self.groups if len(g) > 1]
        if mixed or not ok.all():
            lines.append("")
            lines.append("* not reliable individually: the state mixes with "
                         "others in the field or was matched with a low "
                         "overlap.")
        for g, dg in zip(self.groups, self.group_delta_dipole):
            if len(g) > 1:
                lines.append(
                    "  group " + "+".join(str(n + 1) for n in g)
                    + ": sum of dmu " + "  ".join(
                        f"{'xyz'[i]} {dg[i]:+.4f}" for i in comp))
        return "\n".join(lines)


def clone_options(options):
    """
    Return an independent copy of a pyxtp Options object of dftgwbse.

    Options cannot be deep-copied; the values are transferred through
    their flat dictionary onto fresh defaults.
    """
    from .options import Options

    votcashare = os.environ.get("VOTCASHARE")
    if votcashare is None:
        raise RuntimeError("pyxtp: please set the VOTCASHARE environment "
                           "variable")
    new = Options(f"{votcashare}/xtp/xml/dftgwbse.xml")
    for path, value in options.to_flat_dict().items():
        keys = path.split("/")
        obj = new
        for k in keys[:-1]:
            obj = getattr(obj, k)
        setattr(obj, keys[-1], value)
    return new


@contextlib.contextmanager
def _working_directory(path: Path):
    old = Path.cwd()
    path.mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


class FiniteFieldDipoles:
    """Excited-state dipoles from dftgwbse runs in uniform fields."""

    def __init__(self, atoms, options=None, nthreads: int = 1,
                 field: float = 1e-3, components: str = "xyz",
                 order: int = 2, spin: str = "singlet", nstates: int = 5,
                 extra_states: int = 4, overlap_threshold: float = 0.9,
                 workdir: os.PathLike = "finite_field"):
        """
        Set up the calculation; nothing is run yet.

        Parameters
        ----------
        atoms : ase.Atoms
            The molecule (Angstrom).
        options : pyxtp.options.Options, optional
            Options of ``dftgwbse`` to start from (copied; e.g.
            ``xtp().options`` after setting basis etc.). The field and the
            excited-state settings below are added.
        nthreads : int
            Threads per dftgwbse run.
        field : float
            Step of the field in Hartree/bohr.
        components : str
            Field directions, any combination of "x", "y", "z".
        order : int
            2 (two runs per direction) or 4 (four runs per direction, also
            gives a linearity check).
        spin : str
            "singlet" or "triplet".
        nstates : int
            Number of states to report.
        extra_states : int
            Additional BSE states computed for the matching.
        overlap_threshold : float
            States matched with a lower overlap, or mixing with others, are
            flagged and grouped.
        workdir : os.PathLike
            Directory for the individual runs (kept for inspection).

        """
        if order not in (2, 4):
            raise ValueError("order must be 2 or 4")
        if not components or any(c not in _AXES for c in components):
            raise ValueError("components must be a combination of x, y, z")
        if spin not in ("singlet", "triplet"):
            raise ValueError("spin must be singlet or triplet")
        self.atoms = atoms.copy()
        self.nthreads = nthreads
        self.field = field
        self.components = "".join(dict.fromkeys(components))
        self.order = order
        self.spin = spin
        self.nstates = nstates
        self.extra_states = extra_states
        self.overlap_threshold = overlap_threshold
        self.workdir = Path(workdir).absolute()
        if options is None:
            from .xtp import xtp

            options = xtp().options
        self.options = clone_options(options)
        gwbse = self.options.gwbse
        gwbse.tasks = "gw,singlets" if spin == "singlet" else "gw,triplets"
        gwbse.bse.exctotal = nstates + extra_states
        # converged energies: the field changes Omega by ~1e-4 Hartree
        self.options.dftpackage.convergence_tightness = "verytight"

    def _steps(self) -> List[int]:
        return [-2, -1, 1, 2] if self.order == 4 else [-1, 1]

    def run_point(self, label: str, field_vector: np.ndarray) -> OrbData:
        """One dftgwbse run in the given field (reused if already done)."""
        rundir = self.workdir / label
        orb = rundir / f"{label}.orb"
        if not orb.exists():
            from .xtp import xtp

            opts = clone_options(self.options)
            opts.dftpackage.externalfield = " ".join(
                f"{x:.12g}" for x in field_vector)
            with _working_directory(rundir):
                calc = xtp(options=opts, nthreads=self.nthreads,
                           label=label, directory=".")
                atoms = self.atoms.copy()
                atoms.calc = calc
                calc.calculate_energies(atoms)
        return read_orb(orb, self.spin)

    def run(self) -> FiniteFieldResult:
        """Run all field points and evaluate the derivatives."""
        ref = self.run_point("field_0", np.zeros(3))
        n = min(self.nstates, ref.omega.size)
        nan3 = np.full(3, np.nan)
        dmu = np.full((n, 3), np.nan)
        dalpha = np.full((n, 3), np.nan)
        lin = np.full((n, 3), np.nan) if self.order == 4 else None
        mu_el, alpha0 = nan3.copy(), nan3.copy()
        min_ov = np.ones(n)
        overlaps = []
        runs = {}
        for c in self.components:
            i = _AXES[c]
            for k in self._steps():
                f = np.zeros(3)
                f[i] = k * self.field
                tag = f"{c}{'+' if k > 0 else '-'}{abs(k)}"
                data = self.run_point(f"field_{tag}", f)
                ov = state_overlaps(ref, data)
                cols, vals = match_states(ov, n)
                overlaps.append(ov)
                min_ov = np.minimum(min_ov, vals)
                runs[(c, k)] = (data, cols)

        groups = mixing_groups(overlaps, n, self.overlap_threshold)
        group_dmu = [nan3.copy() for _ in groups]
        for c in self.components:
            i = _AXES[c]
            e = {0: ref.dft_energy}
            om = {0: ref.omega[:n]}
            for k in self._steps():
                data, cols = runs[(c, k)]
                e[k] = data.dft_energy
                om[k] = data.omega[cols]
            mu_el[i] = central_first(e, self.field, self.order)
            alpha0[i] = -central_second(e, self.field, self.order)
            dmu[:, i] = central_first(om, self.field, self.order)
            dalpha[:, i] = -central_second(om, self.field, self.order)
            if lin is not None:
                lin[:, i] = np.abs(
                    central_first(om, self.field, 2)
                    - central_first({k: om[2 * k] for k in (-1, 0, 1)},
                                    2 * self.field, 2))
            for g, gd in zip(groups, group_dmu):
                gd[i] = central_first({k: om[k][g].sum() for k in om},
                                      self.field, self.order)

        com = center_of_mass(ref.elements, ref.positions)
        mu_nuc = (ref.charges[:, None] * (ref.positions - com)).sum(axis=0)
        ground = np.where(np.isnan(mu_el), np.nan, mu_nuc + mu_el)
        return FiniteFieldResult(self.components, self.field, self.order,
                                 self.spin, ref.omega[:n], ground, alpha0,
                                 dmu, dalpha, min_ov, groups, group_dmu, lin,
                                 self.overlap_threshold)


def write_xml(result: FiniteFieldResult, filename: os.PathLike) -> None:
    """Write a FiniteFieldResult as XML."""
    tree = ET.ElementTree(result.to_xml())
    ET.indent(tree, space="  ")
    tree.write(filename, encoding="unicode", xml_declaration=False)
