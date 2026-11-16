"""Finite-difference Hessian from analytic forces.

Every displaced geometry is an independent job in its own directory
(``<workdir>/<name>/``) whose result (forces, energy, positions) is stored
in ``result.npz``. Jobs are claimed with a lock file, so

- an interrupted run resumes where it stopped,
- several processes (on one machine or on cluster nodes sharing the
  file system) can work on the same workdir at the same time, and
- ``run(part=i, nparts=n)`` splits the jobs deterministically, e.g. for
  array jobs.

Any ASE calculator works; ``xtp_ground_state_factory`` sets up pyxtp for
DFT forces with a warm start from the undisplaced geometry.

Central differences (default) need 6N force calculations (+1 reference);
forward differences 3N + 1 but have an O(h) error. The Hessian is
symmetrized; its asymmetry is a measure of the numerical noise.
"""

from __future__ import annotations

import contextlib
import dataclasses
import json
import os
import socket
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Union

import numpy as np
from ase import Atoms

from .modes import NormalModes

Pathlike = Union[str, Path]
_AXES = "xyz"


@dataclasses.dataclass(frozen=True)
class Displacement:
    """One displaced geometry; atom -1 is the reference geometry."""

    name: str
    atom: int
    axis: int
    sign: int

    @classmethod
    def reference(cls) -> "Displacement":
        return cls("eq", -1, -1, 0)


CalculatorFactory = Callable[[Path, Displacement], object]


@contextlib.contextmanager
def _working_directory(path: Path):
    old = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


class FiniteDifferenceHessian:
    """Finite-difference Hessian of an ASE calculator's forces.

    Args:
        atoms: reference geometry (normally an optimized minimum)
        calculator_factory: f(directory, displacement) -> ASE calculator;
            called once per job, with the job directory as the working
            directory during the calculation
        workdir: directory of the jobs and the metadata (hessian.json)
        step: displacement in Angstrom (default 0.005)
        scheme: "central" or "forward"
        indices: atoms to displace (default all); with a subset, the
            Hessian is the partial Hessian of these atoms (the others are
            infinitely heavy)
    """

    def __init__(self, atoms: Atoms,
                 calculator_factory: Optional[CalculatorFactory] = None,
                 workdir: Pathlike = "hessian", step: float = 0.005,
                 scheme: str = "central",
                 indices: Optional[Sequence[int]] = None):
        if scheme not in ("central", "forward"):
            raise ValueError("scheme must be 'central' or 'forward'")
        self.atoms = atoms.copy()
        self.atoms.calc = None
        self.factory = calculator_factory
        self.workdir = Path(workdir).absolute()
        self.step = float(step)
        self.scheme = scheme
        self.indices = list(range(len(atoms))) if indices is None \
            else [int(i) for i in indices]
        self._check_metadata()

    # ------------------------------------------------------------------
    # bookkeeping
    # ------------------------------------------------------------------
    def _metadata(self) -> Dict:
        return {"symbols": self.atoms.get_chemical_symbols(),
                "positions": self.atoms.get_positions().tolist(),
                "step": self.step, "scheme": self.scheme,
                "indices": self.indices}

    def _check_metadata(self) -> None:
        meta_file = self.workdir / "hessian.json"
        meta = self._metadata()
        if meta_file.exists():
            old = json.loads(meta_file.read_text())
            same = (old["symbols"] == meta["symbols"]
                    and np.allclose(old["positions"], meta["positions"],
                                    atol=1e-8)
                    and old["step"] == meta["step"]
                    and old["scheme"] == meta["scheme"]
                    and old["indices"] == meta["indices"])
            if not same:
                raise ValueError(f"{self.workdir} holds a different "
                                 "Hessian (geometry, step, scheme or atoms)"
                                 "; use another workdir")
        else:
            self.workdir.mkdir(parents=True, exist_ok=True)
            tmp = meta_file.with_suffix(f".tmp{os.getpid()}")
            tmp.write_text(json.dumps(meta, indent=1))
            os.replace(tmp, meta_file)

    @classmethod
    def open(cls, workdir: Pathlike,
             calculator_factory: Optional[CalculatorFactory] = None
             ) -> "FiniteDifferenceHessian":
        """Reopen an existing workdir (e.g. to collect or to continue)."""
        meta = json.loads((Path(workdir) / "hessian.json").read_text())
        atoms = Atoms(meta["symbols"], positions=meta["positions"])
        return cls(atoms, calculator_factory, workdir, meta["step"],
                   meta["scheme"], meta["indices"])

    def displacements(self) -> List[Displacement]:
        """All jobs, the reference geometry first."""
        signs = (1, -1) if self.scheme == "central" else (1,)
        jobs = [Displacement.reference()]
        for a in self.indices:
            for k in range(3):
                for s in signs:
                    jobs.append(Displacement(
                        f"a{a:04d}{_AXES[k]}{'+' if s > 0 else '-'}",
                        a, k, s))
        return jobs

    def displaced_atoms(self, d: Displacement) -> Atoms:
        """The geometry of a job."""
        atoms = self.atoms.copy()
        if d.atom >= 0:
            pos = atoms.get_positions()
            pos[d.atom, d.axis] += d.sign * self.step
            atoms.set_positions(pos)
        return atoms

    def _result_file(self, d: Displacement) -> Path:
        return self.workdir / d.name / "result.npz"

    def is_done(self, d: Displacement) -> bool:
        return self._result_file(d).exists()

    def status(self) -> Dict[str, List[str]]:
        """Names of done, running (locked) and pending jobs."""
        out = {"done": [], "running": [], "pending": []}
        for d in self.displacements():
            if self.is_done(d):
                out["done"].append(d.name)
            elif (self.workdir / d.name / ".lock").exists():
                out["running"].append(d.name)
            else:
                out["pending"].append(d.name)
        return out

    def clear_locks(self) -> int:
        """Remove lock files of unfinished jobs (after a crash)."""
        n = 0
        for d in self.displacements():
            lock = self.workdir / d.name / ".lock"
            if lock.exists() and not self.is_done(d):
                lock.unlink()
                n += 1
        return n

    # ------------------------------------------------------------------
    # running
    # ------------------------------------------------------------------
    def _claim(self, d: Displacement) -> bool:
        jobdir = self.workdir / d.name
        jobdir.mkdir(parents=True, exist_ok=True)
        try:
            fd = os.open(jobdir / ".lock",
                         os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            return False
        with os.fdopen(fd, "w") as f:
            f.write(f"{socket.gethostname()} {os.getpid()} {time.time()}\n")
        if self.is_done(d):  # finished by someone else meanwhile
            (jobdir / ".lock").unlink()
            return False
        return True

    def run_one(self, d: Displacement) -> bool:
        """Compute one job unless done or claimed; True if computed."""
        if self.factory is None:
            raise RuntimeError("no calculator_factory given")
        if self.is_done(d) or not self._claim(d):
            return False
        jobdir = self.workdir / d.name
        try:
            atoms = self.displaced_atoms(d)
            with _working_directory(jobdir):
                atoms.calc = self.factory(jobdir, d)
                forces = atoms.get_forces()
                try:
                    energy = atoms.get_potential_energy()
                except Exception:  # some calculators give forces only
                    energy = np.nan
            tmp = jobdir / f"result.tmp{os.getpid()}.npz"
            np.savez(tmp, forces=forces, energy=energy,
                     positions=atoms.get_positions())
            os.replace(tmp, self._result_file(d))
        finally:
            lock = jobdir / ".lock"
            if lock.exists():
                lock.unlink()
        return True

    def run(self, part: int = 0, nparts: int = 1,
            reference_first: bool = True) -> int:
        """Compute all pending jobs (of one part); returns the number run.

        The reference geometry is computed first (by part 0) so that
        factories can warm-start from it; other parts wait for it if
        reference_first is True.
        """
        jobs = self.displacements()
        ref, rest = jobs[0], jobs[1:]
        count = 0
        if part == 0:
            count += self.run_one(ref)
        elif reference_first:
            while not self.is_done(ref):
                time.sleep(5)
        for d in rest[part::nparts]:
            count += self.run_one(d)
        return count

    # ------------------------------------------------------------------
    # results
    # ------------------------------------------------------------------
    def _load(self, d: Displacement):
        f = self._result_file(d)
        if not f.exists():
            raise RuntimeError(f"job {d.name} is not done; "
                               f"status: {self._count()}")
        with np.load(f) as r:
            return r["forces"], float(r["energy"])

    def _count(self) -> str:
        st = self.status()
        return ", ".join(f"{k} {len(v)}" for k, v in st.items())

    def reference_forces(self) -> np.ndarray:
        """Forces at the reference geometry (eV/Angstrom)."""
        return self._load(Displacement.reference())[0]

    def hessian(self, symmetrize: bool = True) -> np.ndarray:
        """Hessian (3n x 3n, eV/Angstrom^2) over the displaced atoms."""
        idx = np.array(self.indices)
        rows = np.ravel([[3 * a, 3 * a + 1, 3 * a + 2] for a in idx])
        h = np.zeros((rows.size, rows.size))
        jobs = {(d.atom, d.axis, d.sign): d for d in self.displacements()}
        f0 = None
        if self.scheme == "forward":
            f0 = self.reference_forces().ravel()[rows]
        for col, (a, k) in enumerate((a, k) for a in idx for k in range(3)):
            fp = self._load(jobs[(a, k, 1)])[0].ravel()[rows]
            if self.scheme == "central":
                fm = self._load(jobs[(a, k, -1)])[0].ravel()[rows]
                h[:, col] = -(fp - fm) / (2 * self.step)
            else:
                h[:, col] = -(fp - f0) / self.step
        self.asymmetry = float(np.abs(h - h.T).max())
        if symmetrize:
            h = 0.5 * (h + h.T)
        return h

    def normal_modes(self) -> NormalModes:
        """Normal modes of the (partial) Hessian."""
        h = self.hessian()
        sub = self.atoms[self.indices]
        return NormalModes.from_hessian(sub, h)


# ----------------------------------------------------------------------
# pyxtp
# ----------------------------------------------------------------------
def xtp_ground_state_factory(settings: Optional[Dict[str, str]] = None,
                             nthreads: int = 1,
                             warm_start: bool = True,
                             options=None) -> CalculatorFactory:
    """Factory for pyxtp ground-state DFT forces.

    Runs only the DFT part of dftgwbse (tasks input,dft,parse) with
    analytic forces and verytight SCF convergence. With warm_start, the
    displaced geometries start from the converged MOs of the reference
    job (moguess + initial_guess orbfile).

    Args:
        settings: pyxtp options, e.g.
            {"dftpackage/functional": "XC_HYB_GGA_XC_PBEH",
            "dftpackage/basisset": "def2-svp",
            "dftpackage/auxbasisset": "aux-def2-svp"}; they override the
            defaults set here
        nthreads: threads per calculation
        warm_start: start displaced SCFs from the reference MOs
        options: optional pyxtp Options object to start from
    """
    from ..xtp import xtp  # needs the compiled bindings

    def factory(jobdir: Path, d: Displacement):
        # the runner calls the calculator inside jobdir; pyxtp's own
        # 'directory' would change its file names (ASE label), so it stays
        # at the default
        calc = xtp(nthreads=nthreads)
        if options is not None:
            # copy the values (calc.set changes calc.options in place, and
            # deepcopy of Options does not work)
            calc.set(**options.to_flat_dict())
        base = {"tasks": "input,dft,parse",
                "dftpackage/xtpdft/compute_forces": "true",
                "dftpackage/convergence_tightness": "verytight"}
        if warm_start and d is not None and d.atom >= 0:
            orbs = sorted((jobdir.parent / "eq").glob("*.orb"))
            if orbs:
                base["moguess"] = orbs[0].as_posix()
                base["dftpackage/initial_guess"] = "orbfile"
        base.update(settings or {})
        calc.set(**base)
        return calc

    return factory
