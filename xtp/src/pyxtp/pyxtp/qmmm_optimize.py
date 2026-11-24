"""Geometry optimization of the QM region of a QM/MM (qmmm) job.

Every energy and force evaluation is one ``xtp_parallel -e qmmm`` job in its
own step directory, with

* the QM coordinates of the step in the frame of the job
  (``qmregion.geometry``, instead of the segment mapping),
* the final orbitals of the previous step as the first SCF guess
  (``qmregion.guess``),
* the induced moments of the previous step as the start of the polar
  regions (``polarregion.guess``, the last checkpoint of that job),
* ``qmregion.compute_forces`` on, so that the job output holds the energy
  (``E_tot``), the forces on the QM atoms and the coordinates they belong to.

The environment is frozen at its topology: the region partition and the
Ewald background stay as set up from the segment positions, and the polar
regions are re-converged in every job. For excited states set
``qmregion.relaxed_density`` so that the energy is stationary in the induced
moments and the forces are its derivatives.

The first evaluation runs the job as given (mapped geometry) and provides the
start geometry in the job frame (``QMMMCalculator.start_atoms``). Any ASE
optimizer can then drive the calculator; ``optimize`` wraps BFGS, FIRE or
Sella.

Terms added in python (``terms``): the qmmm energy of a QM region in
point-multipole/polarizable embedding has no short-range repulsion between
QM and environment atoms, and xtp's DFT has no dispersion. In an
optimization a QM atom can then be pulled onto an environment charge (e.g.
the hydroxyl H of a solute onto the O of a water). ``LennardJones`` adds
Lennard-Jones interactions between QM and environment atoms (as the
QM-MM van der Waals term of standard QM/MM, parameters typically those of
the force field of the MD), ``D3Dispersion`` (simple-dftd3, the python
package ``dftd3``) Grimme's D3 within the QM region and, optionally, with
the environment. The environment atoms are taken from the checkpoint of
the first job (sites with an element name) and stay frozen; each term
keeps the set of environment atoms within its cutoff at the start, so
energy and forces stay smooth. For D3 with the environment the energy is
E_D3(QM + environment) - E_D3(environment), so the forces on the QM atoms
are exact (the coordination numbers of the environment atoms, and the
three-body term if switched on, follow the QM atoms). Use either the
dispersion of LJ or that of D3 between QM and environment, not both:
with LJ, D3 within the QM region only.

Units: ASE units (eV, Angstrom) towards the optimizer; the job output gives
E_tot in eV and the forces in Hartree/Bohr.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from typing import List, Optional, Sequence

import numpy as np
from ase import Atoms, units
from ase.calculators.calculator import Calculator, all_changes

# entries of the setup directory that are not linked into the step directories
_NOT_LINKED = {"QMMM"}


def environment_from_checkpoint(checkpoint: str):
    """Atoms of the polar and static regions of a qmmm checkpoint: element
    symbols and positions (Angstrom, job frame), and the QM atoms (symbols,
    positions) of the same file for a frame check. Sites whose name is not
    an element (virtual sites) are skipped."""
    import h5py
    from ase.data import atomic_numbers

    env_symbols: List[str] = []
    env_positions: List[List[float]] = []
    qm_symbols: List[str] = []
    qm_positions: List[List[float]] = []

    def text(v):
        return v.decode() if isinstance(v, bytes) else str(v)

    with h5py.File(checkpoint, "r") as h5:
        for name, region in h5.items():
            rtype = text(np.ravel(region.attrs.get("type", [""]))[0])
            if rtype == "qmregion":
                atoms = region["QMdata/qmmolecule/qmatoms"][()].ravel()
                for a in atoms:
                    qm_symbols.append(text(a["element"]))
                    qm_positions.append([a["posX"], a["posY"], a["posZ"]])
                continue
            if rtype not in ("polarregion", "staticregion"):
                continue
            dataset = "polarsites" if rtype == "polarregion" else "staticsites"
            for seg in region["segments"].values():
                for site in seg[dataset][()].ravel():
                    el = text(site["type"])
                    if atomic_numbers.get(el, 0) == 0:  # also ASE's dummy X
                        continue
                    env_symbols.append(el)
                    env_positions.append([site["posX"], site["posY"], site["posZ"]])
    bohr = units.Bohr
    return (env_symbols, np.array(env_positions).reshape(-1, 3) * bohr,
            qm_symbols, np.array(qm_positions).reshape(-1, 3) * bohr)


class D3Dispersion:
    """Grimme's D3 for the QM region of a qmmm job (simple-dftd3).

    Parameters
    ----------
    method : functional the damping parameters belong to (e.g. "pbe0").
    damping : "d3bj" (rational, default), "d3zero", "d3bjm" or "d3zerom".
    atm : include the three-body Axilrod-Teller-Muto term (default off,
        as in dftd3.ase and the usual D3(BJ)).
    environment : include the atoms of the polar and static regions (frozen).
    cutoff : environment atoms farther than this from every QM atom are left
        out (Angstrom; None: all).
    """

    def __init__(self, method: str = "pbe0", damping: str = "d3bj",
                 atm: bool = False, environment: bool = True,
                 cutoff: Optional[float] = 15.0):
        try:
            from dftd3 import interface
        except ImportError as err:
            raise ImportError("D3Dispersion needs the python package dftd3 "
                              "(simple-dftd3)") from err
        params = {"d3bj": interface.RationalDampingParam,
                  "d3zero": interface.ZeroDampingParam,
                  "d3bjm": interface.ModifiedRationalDampingParam,
                  "d3zerom": interface.ModifiedZeroDampingParam}
        if damping not in params:
            raise ValueError(f"unknown D3 damping {damping}")
        self.name = "d3"
        self._interface = interface
        self.param = params[damping](method=method, atm=atm)
        self.method = method
        self.damping = damping
        self.environment = environment
        self.cutoff = cutoff
        self.env_symbols: List[str] = []
        self.env_positions = np.zeros((0, 3))
        self._env_energy: Optional[float] = None

    def set_environment(self, symbols: Sequence[str], positions: np.ndarray,
                        qm_positions: Optional[np.ndarray] = None) -> None:
        """Environment atoms (Angstrom); with qm_positions and a cutoff only
        those within the cutoff of a QM atom are kept."""
        positions = np.asarray(positions, dtype=float).reshape(-1, 3)
        keep = np.ones(len(symbols), dtype=bool)
        if self.cutoff is not None and qm_positions is not None and len(symbols):
            d = np.linalg.norm(positions[:, None, :] -
                               np.asarray(qm_positions)[None, :, :], axis=2)
            keep = d.min(axis=1) <= self.cutoff
        self.env_symbols = [s for s, k in zip(symbols, keep) if k]
        self.env_positions = positions[keep]
        self._env_energy = None

    def _model(self, symbols, positions):
        from ase.data import atomic_numbers
        numbers = np.array([atomic_numbers[s] for s in symbols])
        return self._interface.DispersionModel(numbers, positions / units.Bohr)

    def env_energy(self) -> float:
        """D3 energy of the environment alone (eV), a constant."""
        if self._env_energy is None:
            if len(self.env_symbols) == 0:
                self._env_energy = 0.0
            else:
                res = self._model(self.env_symbols, self.env_positions).get_dispersion(
                    self.param, grad=False)
                self._env_energy = res["energy"] * units.Hartree
        return self._env_energy

    def __call__(self, symbols: Sequence[str], positions: np.ndarray):
        """D3 energy (eV) of the QM atoms and their interaction with the
        environment, and the forces on the QM atoms (eV/Angstrom)."""
        nqm = len(symbols)
        syms = list(symbols)
        pos = np.asarray(positions, dtype=float).reshape(-1, 3)
        if self.environment and len(self.env_symbols):
            syms = syms + self.env_symbols
            pos = np.vstack([pos, self.env_positions])
        res = self._model(syms, pos).get_dispersion(self.param, grad=True)
        energy = res["energy"] * units.Hartree
        if self.environment and len(self.env_symbols):
            energy -= self.env_energy()
        forces = -res["gradient"][:nqm] * units.Hartree / units.Bohr
        return energy, forces


class QMMMJobError(RuntimeError):
    """An xtp_parallel qmmm job did not complete."""


def _vector(text: str) -> List[float]:
    return [float(v) for v in text.split()]


def read_job_output(jobfile: str):
    """Energy (eV), forces (Hartree/Bohr), elements, positions (Angstrom),
    initial and final state of the single job in ``jobfile``."""
    root = ET.parse(jobfile).getroot()
    job = root.find("job")
    if job is None:
        raise QMMMJobError(f"{jobfile}: no job")
    status = job.findtext("status", "").strip()
    if status != "COMPLETE":
        raise QMMMJobError(f"{jobfile}: job status {status}")
    output = job.find("output")
    region = None
    for r in output.iter("region"):
        if r.get("type") == "qmregion":
            region = r
    if region is None:
        raise QMMMJobError(f"{jobfile}: no qmregion in the job output")
    forces = region.find("forces")
    geometry = region.find("geometry")
    if forces is None or geometry is None:
        raise QMMMJobError(
            f"{jobfile}: no forces/geometry in the job output (needs "
            "qmregion.compute_forces and an xtp version with the qmregion "
            "geometry output)")
    if geometry.get("unit", "Angstrom") != "Angstrom":
        raise QMMMJobError(f"{jobfile}: unexpected geometry unit")
    if forces.get("unit", "Hartree/Bohr") != "Hartree/Bohr":
        raise QMMMJobError(f"{jobfile}: unexpected force unit")
    elements = [a.get("element") for a in geometry.findall("atom")]
    positions = np.array([_vector(a.text) for a in geometry.findall("atom")])
    f = np.array([_vector(a.text) for a in forces.findall("atom")])
    if f.shape != positions.shape:
        raise QMMMJobError(f"{jobfile}: forces and geometry do not match")
    return {
        "energy": float(output.findtext("E_tot")),
        "forces": f,
        "elements": elements,
        "positions": positions,
        "initial_state": region.findtext("Initial_State", "").strip(),
        "final_state": region.findtext("Final_State", "").strip(),
        "forces_state": forces.get("state", ""),
    }


def _find(directory: str, name: str) -> Optional[str]:
    for dirpath, _, files in os.walk(directory):
        if name in files:
            return os.path.join(dirpath, name)
    return None


def _last_checkpoint(directory: str) -> Optional[str]:
    """checkpoint_iter_N.hdf5 with the largest N below directory."""
    best, best_n = None, -1
    for dirpath, _, files in os.walk(directory):
        for name in files:
            m = re.fullmatch(r"checkpoint_iter_(\d+)\.hdf5", name)
            if m and int(m.group(1)) > best_n:
                best, best_n = os.path.join(dirpath, name), int(m.group(1))
    return best


class LennardJones:
    """Lennard-Jones interaction between the QM atoms and the (frozen)
    environment atoms, 4 eps ((s/r)^12 - (s/r)^6), Lorentz-Berthelot mixing.

    Parameters
    ----------
    qm : per QM atom, either a dict element -> (sigma, epsilon) or a list
        of (sigma, epsilon) in the order of the QM atoms.
    environment : dict element -> (sigma, epsilon) for the environment
        atoms; elements that are missing get no LJ (warning once).
    cutoff : environment atoms farther than this from every QM atom at the
        start are left out (Angstrom; None: all).
    Units: sigma in Angstrom, epsilon in kJ/mol.
    """

    def __init__(self, qm, environment: dict, cutoff: Optional[float] = 12.0):
        self.name = "lj"
        self.environment = True
        self.qm = qm
        self.env_params = {k: tuple(v) for k, v in environment.items()}
        self.cutoff = cutoff
        self.env_symbols: List[str] = []
        self.env_positions = np.zeros((0, 3))
        self._sig = np.zeros(0)
        self._eps = np.zeros(0)
        self.missing: List[str] = []

    def set_environment(self, symbols: Sequence[str], positions: np.ndarray,
                        qm_positions: Optional[np.ndarray] = None) -> None:
        positions = np.asarray(positions, dtype=float).reshape(-1, 3)
        keep = np.array([s in self.env_params for s in symbols], dtype=bool)
        self.missing = sorted({s for s in symbols if s not in self.env_params})
        if self.cutoff is not None and qm_positions is not None and len(symbols):
            d = np.linalg.norm(positions[:, None, :] -
                               np.asarray(qm_positions)[None, :, :], axis=2)
            keep &= d.min(axis=1) <= self.cutoff
        self.env_symbols = [s for s, k in zip(symbols, keep) if k]
        self.env_positions = positions[keep]
        self._sig = np.array([self.env_params[s][0] for s in self.env_symbols])
        self._eps = np.array([self.env_params[s][1] for s in self.env_symbols])

    def _qm_params(self, symbols):
        if isinstance(self.qm, dict):
            missing = sorted({s for s in symbols if s not in self.qm})
            if missing:
                raise ValueError(f"LennardJones: no parameters for QM "
                                 f"element(s) {missing}")
            return (np.array([self.qm[s][0] for s in symbols]),
                    np.array([self.qm[s][1] for s in symbols]))
        if len(self.qm) != len(symbols):
            raise ValueError(f"LennardJones: {len(self.qm)} parameter pairs for "
                             f"{len(symbols)} QM atoms")
        p = np.asarray(self.qm, dtype=float)
        return p[:, 0], p[:, 1]

    def __call__(self, symbols: Sequence[str], positions: np.ndarray):
        """Energy (eV) and forces on the QM atoms (eV/Angstrom)."""
        pos = np.asarray(positions, dtype=float).reshape(-1, 3)
        if len(self.env_symbols) == 0:
            return 0.0, np.zeros_like(pos)
        sq, eq = self._qm_params(symbols)
        sig = 0.5 * (sq[:, None] + self._sig[None, :])
        eps = np.sqrt(eq[:, None] * self._eps[None, :])
        rvec = pos[:, None, :] - self.env_positions[None, :, :]
        r = np.linalg.norm(rvec, axis=2)
        sr6 = (sig / r) ** 6
        kj = units.kJ / units.mol  # kJ/mol in eV
        energy = (4.0 * eps * (sr6 * sr6 - sr6)).sum() * kj
        # -dE/dr / r times rvec
        coef = 4.0 * eps * (12.0 * sr6 * sr6 - 6.0 * sr6) / (r * r)
        forces = (coef[:, :, None] * rvec).sum(axis=1) * kj
        return energy, forces


class QMMMCalculator(Calculator):
    """ASE calculator: energy and forces of the QM region of a qmmm job.

    Parameters
    ----------
    setup : directory with the qmmm setup (state file, mapping, MD/QC/MP
        files, ...); its entries are linked into every step directory.
    options : qmmm options file (``<options><qmmm>...``), relative to setup.
    jobfile : job file with the job to run, relative to setup.
    job_id : id of the job in ``jobfile`` (default: the first job).
    statefile : state file passed to ``xtp_parallel -f``.
    workdir : directory for the step directories.
    executable : ``xtp_parallel``.
    threads : threads per job (``-t``).
    terms : terms added to energy and forces (``LennardJones``,
        ``D3Dispersion``), each with a distinct ``name``.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, setup: str = ".", options: str = "qmmm.xml",
                 jobfile: Optional[str] = None, job_id: Optional[str] = None,
                 statefile: str = "state.hdf5", workdir: str = "qmmm_opt",
                 executable: str = "xtp_parallel", threads: int = 1,
                 terms: Sequence = (), **kwargs):
        super().__init__(**kwargs)
        self.setup = os.path.abspath(setup)
        self.options_file = os.path.join(self.setup, options)
        self.options_tree = ET.parse(self.options_file)
        qmmm = self.options_tree.getroot().find("qmmm")
        if qmmm is None:
            raise ValueError(f"{self.options_file}: expected <options><qmmm>")
        if jobfile is None:
            jobfile = qmmm.findtext("job_file", "").strip()
            if not jobfile:
                raise ValueError("no job file given or set in the options")
        self.jobfile = os.path.join(self.setup, jobfile)
        self.job = self._select_job(job_id)
        self.statefile = statefile
        self.workdir = os.path.abspath(workdir)
        self.executable = executable
        self.threads = threads
        self.nsteps = 0
        self.history: List[dict] = []
        self._last_orb: Optional[str] = None
        self._last_checkpoint: Optional[str] = None
        self._start = None
        self.terms = list(terms)
        names = [t.name for t in self.terms]
        if len(set(names)) != len(names):
            raise ValueError(f"terms with the same name: {names}")
        self._terms_ready = False
        os.makedirs(self.workdir, exist_ok=True)

    # ---- setup -------------------------------------------------------------
    def _select_job(self, job_id):
        jobs = ET.parse(self.jobfile).getroot().findall("job")
        if not jobs:
            raise ValueError(f"{self.jobfile}: no jobs")
        if job_id is None:
            job = jobs[0]
        else:
            matches = [j for j in jobs if j.findtext("id", "").strip() == str(job_id)]
            if len(matches) != 1:
                raise ValueError(f"{self.jobfile}: no unique job with id {job_id}")
            job = matches[0]
        job = ET.fromstring(ET.tostring(job))
        for tag in ("output", "host", "time"):
            for child in job.findall(tag):
                job.remove(child)
        status = job.find("status")
        if status is None:
            status = ET.SubElement(job, "status")
        status.text = "AVAILABLE"
        return job

    def _write_options(self, stepdir: str, geometry: bool) -> None:
        tree = ET.ElementTree(ET.fromstring(ET.tostring(self.options_tree.getroot())))
        qmmm = tree.getroot().find("qmmm")
        jf = qmmm.find("job_file")
        if jf is None:
            jf = ET.SubElement(qmmm, "job_file")
        jf.text = "jobs.xml"
        qmregion = qmmm.find("regions/qmregion")
        if qmregion is None:
            raise ValueError(f"{self.options_file}: no regions/qmregion")
        for tag in ("compute_forces", "geometry", "guess"):
            for child in qmregion.findall(tag):
                qmregion.remove(child)
        ET.SubElement(qmregion, "compute_forces").text = "true"
        if geometry:
            ET.SubElement(qmregion, "geometry").text = os.path.join(stepdir, "geometry.xyz")
        if self._last_orb is not None:
            ET.SubElement(qmregion, "guess").text = self._last_orb
        for polar in qmmm.findall("regions/polarregion"):
            for child in polar.findall("guess"):
                polar.remove(child)
            if self._last_checkpoint is not None:
                ET.SubElement(polar, "guess").text = self._last_checkpoint
        tree.write(os.path.join(stepdir, "qmmm.xml"))

    def _prepare_step(self, stepdir: str) -> None:
        if os.path.exists(stepdir):
            shutil.rmtree(stepdir)
        os.makedirs(stepdir)
        skip = {os.path.basename(self.options_file), os.path.basename(self.jobfile)}
        for entry in os.listdir(self.setup):
            if entry in skip or entry in _NOT_LINKED:
                continue
            src = os.path.join(self.setup, entry)
            if os.path.abspath(src) == self.workdir:
                continue
            os.symlink(src, os.path.join(stepdir, entry))
        root = ET.Element("jobs")
        root.append(ET.fromstring(ET.tostring(self.job)))
        ET.ElementTree(root).write(os.path.join(stepdir, "jobs.xml"))

    # ---- running -----------------------------------------------------------
    def run_job(self, positions: Optional[np.ndarray] = None,
                elements: Optional[Sequence[str]] = None) -> dict:
        """One qmmm job; positions (Angstrom, job frame) None: mapped
        geometry. Returns the parsed job output."""
        stepdir = os.path.join(self.workdir, f"step_{self.nsteps:03d}")
        self._prepare_step(stepdir)
        if positions is not None:
            with open(os.path.join(stepdir, "geometry.xyz"), "w") as f:
                f.write(f"{len(elements)}\nqmmm_optimize step {self.nsteps}\n")
                for el, r in zip(elements, positions):
                    f.write(f"{el} {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}\n")
        self._write_options(stepdir, positions is not None)
        cmd = [self.executable, "-e", "qmmm", "-o", "qmmm.xml", "-f",
               self.statefile, "-x", "1", "-t", str(self.threads), "-j", "run"]
        with open(os.path.join(stepdir, "run.log"), "w") as log:
            proc = subprocess.run(cmd, cwd=stepdir, stdout=log,
                                  stderr=subprocess.STDOUT)
        try:
            if proc.returncode != 0:
                raise QMMMJobError(f"xtp_parallel exited with {proc.returncode}")
            result = read_job_output(os.path.join(stepdir, "jobs.xml"))
        except QMMMJobError as err:
            raise QMMMJobError(f"step {self.nsteps} in {stepdir}: {err} "
                               "(see run.log there)") from err
        orb = _find(os.path.join(stepdir, "QMMM"), "qmregion_final.orb")
        if orb is not None:
            self._last_orb = orb
        checkpoint = _last_checkpoint(os.path.join(stepdir, "QMMM"))
        if checkpoint is not None:
            self._last_checkpoint = checkpoint
        if elements is not None and list(elements) != result["elements"]:
            raise QMMMJobError(f"step {self.nsteps}: the elements of the job "
                               "differ from those of the atoms")
        if positions is not None:
            dev = np.abs(result["positions"] - positions).max()
            if dev > 1e-6:
                raise QMMMJobError(f"step {self.nsteps}: the job ran at a "
                                   f"geometry {dev} Angstrom off the given one")
        if result["initial_state"] and result["final_state"] != result["initial_state"]:
            result["warning"] = (f"tracked state changed from "
                                 f"{result['initial_state']} to {result['final_state']}")
        if self.terms:
            self._add_terms(result, checkpoint)
        result["step"] = self.nsteps
        result["stepdir"] = stepdir
        self.history.append(result)
        self.nsteps += 1
        return result

    def _add_terms(self, result: dict, checkpoint: Optional[str]) -> None:
        """The python terms on top of the job's energy and forces; the
        environment comes from the first job's checkpoint (it is frozen)."""
        if not self._terms_ready and any(t.environment for t in self.terms):
            if checkpoint is None:
                raise QMMMJobError(f"step {self.nsteps}: no checkpoint to take "
                                   "the environment atoms from")
            env_s, env_p, qm_s, qm_p = environment_from_checkpoint(checkpoint)
            if qm_s != result["elements"] or (
                    len(qm_p) and np.abs(qm_p - result["positions"]).max() > 1e-4):
                raise QMMMJobError(f"step {self.nsteps}: the QM atoms of the "
                                   "checkpoint differ from the job output; "
                                   "cannot place the environment")
            for t in self.terms:
                if t.environment:
                    t.set_environment(env_s, env_p, result["positions"])
        self._terms_ready = True
        result["energy_qmmm"] = result["energy"]
        for t in self.terms:
            e, f = t(result["elements"], result["positions"])
            result[f"energy_{t.name}"] = e
            result["energy"] = result["energy"] + e
            result["forces"] = result["forces"] + f / (units.Hartree / units.Bohr)

    def start_atoms(self) -> Atoms:
        """The QM region at the mapped geometry, in the job frame (runs the
        first job; its energy and forces are reused for these positions)."""
        if self._start is None:
            self._start = self.run_job()
        atoms = Atoms(self._start["elements"], positions=self._start["positions"])
        atoms.calc = self
        return atoms

    def calculate(self, atoms=None, properties=("energy", "forces"),
                  system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        pos = self.atoms.get_positions()
        res = None
        if (self._start is not None and self.nsteps == 1
                and self._start["elements"] == self.atoms.get_chemical_symbols()
                and np.allclose(self._start["positions"], pos, atol=1e-9, rtol=0)):
            res = self._start
        if res is None:
            res = self.run_job(pos, self.atoms.get_chemical_symbols())
        self.results["energy"] = res["energy"]
        self.results["forces"] = res["forces"] * units.Hartree / units.Bohr


def optimize(calc: QMMMCalculator, fmax: float = 0.05, steps: int = 100,
             optimizer: str = "BFGS", fix: Sequence[int] = (),
             trajectory: Optional[str] = None, logfile: Optional[str] = "-",
             atoms: Optional[Atoms] = None, maxstep: Optional[float] = None,
             rigid_body: bool = True):
    """Optimize the QM region with an ASE optimizer (BFGS, FIRE or Sella).

    fmax in eV/Angstrom; fix: indices of QM atoms kept in place; maxstep:
    largest atomic step in Angstrom (BFGS, FIRE). rigid_body False (Sella
    only): the position and orientation of the QM region are kept, only its
    internal coordinates relax; fmax then refers to the forces without
    their net force and torque (a frozen cage that does not fit the QM
    model can otherwise drag the molecule along for many steps). Returns
    the optimized atoms (job frame) and whether fmax was reached."""
    from ase.constraints import FixAtoms

    if atoms is None:
        atoms = calc.start_atoms()
    atoms.calc = calc
    if fix:
        atoms.set_constraint(FixAtoms(indices=list(fix)))
    name = optimizer.lower()
    if not rigid_body and name != "sella":
        raise ValueError("rigid_body=False needs the Sella optimizer")
    kwargs = {} if maxstep is None else {"maxstep": maxstep}
    if name == "bfgs":
        from ase.optimize import BFGS as Opt
        opt = Opt(atoms, trajectory=trajectory, logfile=logfile, **kwargs)
    elif name == "fire":
        from ase.optimize import FIRE as Opt
        opt = Opt(atoms, trajectory=trajectory, logfile=logfile, **kwargs)
    elif name == "sella":
        # The QM region sits in a fixed environment: its energy is not
        # invariant under translations and rotations, which Sella projects
        # out by default (as for a free molecule). Its "converged" geometry
        # would then keep a net force and torque. Internal coordinates with
        # explicit translations and rotations of the whole region instead.
        from sella import Internals, Sella
        ints = Internals(atoms)
        ints.find_all_bonds()
        ints.find_all_angles()
        ints.find_all_dihedrals()
        if rigid_body:
            ints.add_translation(tuple(range(len(atoms))))
            ints.add_rotation(tuple(range(len(atoms))))
            opt = Sella(atoms, order=0, internal=ints, trajectory=trajectory,
                        logfile=logfile)
        else:
            # Sella's free-molecule default: translations and rotations
            # projected out
            opt = Sella(atoms, order=0, internal=True, trajectory=trajectory,
                        logfile=logfile)
    else:
        raise ValueError(f"unknown optimizer {optimizer}")
    converged = bool(opt.run(fmax=fmax, steps=steps))
    return atoms, converged
