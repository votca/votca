"""Tests of pyxtp.vibronic (no VOTCA calculations needed)."""
import math
import os
import unittest

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

try:
    from pyxtp.vibronic import (FiniteDifferenceHessian, NormalModes,
                                VibronicCouplings, couplings_from_geometry,
                                couplings_from_gradient, ensemble_spectrum,
                                read_gap_gradient, read_orca_engrad,
                                read_orca_hessian, vibronic_lineshape)
    from pyxtp.vibronic.lineshape import PROCESSES
    from pyxtp.vibronic.modes import KB_EV
except ImportError:
    raise unittest.SkipTest('pyxtp not available')


# ----------------------------------------------------------------------
# a molecule with an analytic Hessian: harmonic springs between all pairs
# ----------------------------------------------------------------------
class Springs(Calculator):
    """E = sum_{i<j} k_ij/2 (r_ij - r0_ij)^2 (eV, Angstrom)."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, k, r0, **kwargs):
        super().__init__(**kwargs)
        self.k = k
        self.r0 = r0

    def calculate(self, atoms=None, properties=None,
                  system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        p = self.atoms.get_positions()
        n = len(p)
        e = 0.0
        f = np.zeros_like(p)
        for i in range(n):
            for j in range(i + 1, n):
                d = p[i] - p[j]
                r = np.linalg.norm(d)
                e += 0.5 * self.k[i, j] * (r - self.r0[i, j]) ** 2
                fij = -self.k[i, j] * (r - self.r0[i, j]) * d / r
                f[i] += fij
                f[j] -= fij
        self.results = {"energy": e, "forces": f}


def spring_hessian(p, k, r0):
    n = len(p)
    h = np.zeros((3 * n, 3 * n))
    for i in range(n):
        for j in range(i + 1, n):
            d = p[i] - p[j]
            r = np.linalg.norm(d)
            u = d / r
            b = k[i, j] * (np.outer(u, u) + (1 - r0[i, j] / r)
                           * (np.eye(3) - np.outer(u, u)))
            for a, c, s in ((i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)):
                h[3 * a:3 * a + 3, 3 * c:3 * c + 3] += s * b
    return h


def water_like(seed=0):
    rng = np.random.default_rng(seed)
    atoms = Atoms("OH2CN", positions=[[0, 0, 0], [0.96, 0, 0],
                                      [-0.24, 0.93, 0], [0.3, -0.5, 1.2],
                                      [1.1, -0.7, 1.9]])
    p = atoms.get_positions()
    n = len(atoms)
    k = rng.uniform(5.0, 40.0, (n, n))
    k = 0.5 * (k + k.T)
    r0 = np.linalg.norm(p[:, None] - p[None], axis=-1)
    return atoms, k, r0


def random_rotation(rng):
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q * np.sign(np.linalg.det(q))


# ----------------------------------------------------------------------
# normal modes
# ----------------------------------------------------------------------
def test_modes_count_and_rigid_body_projection():
    atoms, k, r0 = water_like()
    h = spring_hessian(atoms.get_positions(), k, r0)
    m = NormalModes.from_hessian(atoms, h)
    assert m.nmodes == 3 * len(atoms) - 6
    assert (m.eigenvalues > 0).all()
    # rigidly moved molecule: same frequencies
    rng = np.random.default_rng(3)
    rot = random_rotation(rng)
    moved = atoms.copy()
    moved.set_positions(atoms.get_positions() @ rot.T + [1.0, -2.0, 0.5])
    big = np.kron(np.eye(len(atoms)), rot)
    m2 = NormalModes.from_hessian(moved, big @ h @ big.T)
    np.testing.assert_allclose(m2.wavenumbers, m.wavenumbers, rtol=1e-10)


def test_linear_molecule_has_3n_minus_5_modes():
    atoms = Atoms("CO2", positions=[[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]])
    p = atoms.get_positions()
    k = np.full((3, 3), 20.0)
    r0 = np.linalg.norm(p[:, None] - p[None], axis=-1)
    m = NormalModes.from_hessian(atoms, spring_hessian(p, k, r0))
    assert m.nmodes == 4


def test_rigid_body_gradient_has_no_coupling():
    atoms, k, r0 = water_like()
    m = NormalModes.from_hessian(
        atoms, spring_hessian(atoms.get_positions(), k, r0))
    masses = atoms.get_masses()
    # uniform field force ~ mass (pure translation) and a torque
    g = np.outer(masses, [0.1, -0.2, 0.05])
    c = couplings_from_gradient(m, g)
    assert c.total_reorganization < 1e-20
    assert c.net_force > 0.1


def test_gradient_and_geometry_routes_agree_for_harmonic_surface():
    atoms, k, r0 = water_like(1)
    h = spring_hessian(atoms.get_positions(), k, r0)
    m = NormalModes.from_hessian(atoms, h)
    rng = np.random.default_rng(7)
    g = rng.normal(scale=0.3, size=(len(atoms), 3))
    c1 = couplings_from_gradient(m, g)
    # displaced minimum in mass-weighted normal coordinates
    sqm = np.repeat(np.sqrt(m.masses), 3)
    q = m.modes.T @ (g.ravel() / sqm)
    dx = m.modes @ (-q / m.eigenvalues) / sqm
    final = atoms.copy()
    final.set_positions(atoms.get_positions() + dx.reshape(-1, 3))
    rot = random_rotation(rng)
    final.set_positions(final.get_positions() @ rot.T + 3.0)
    c2 = couplings_from_geometry(m, final)
    # the fit is not exactly the Eckart frame: tiny second-order deviation
    np.testing.assert_allclose(c2.reorganization, c1.reorganization,
                               rtol=1e-3, atol=1e-7)


def test_save_load_roundtrip(tmp_path):
    atoms, k, r0 = water_like()
    m = NormalModes.from_hessian(
        atoms, spring_hessian(atoms.get_positions(), k, r0))
    m.save(tmp_path / "m.npz")
    m2 = NormalModes.load(tmp_path / "m.npz")
    np.testing.assert_allclose(m2.wavenumbers, m.wavenumbers)


def test_orca_readers(tmp_path):
    h = np.arange(36, dtype=float).reshape(6, 6)
    h = h + h.T
    lines = ["$orca_hessian_file", "", "$hessian", "6"]
    for c0 in range(0, 6, 5):
        cols = list(range(c0, min(c0 + 5, 6)))
        lines.append("    " + "  ".join(str(c) for c in cols))
        for r in range(6):
            lines.append(f"{r:4d} " + " ".join(f"{h[r, c]:.10E}"
                                               for c in cols))
    lines += ["", "$end"]
    (tmp_path / "x.hess").write_text("\n".join(lines))
    hartree_bohr2 = 27.211386245988 / 0.529177210903**2
    np.testing.assert_allclose(read_orca_hessian(tmp_path / "x.hess"),
                               h * hartree_bohr2, rtol=1e-6)
    eng = "#\n# Number of atoms\n#\n2\n#\n# energy\n#\nn/a\n#\n# grad\n#\n" + \
        "\n".join(str(x) for x in range(6)) + "\n"
    (tmp_path / "x.engrad").write_text(eng)
    g = read_orca_engrad(tmp_path / "x.engrad")
    np.testing.assert_allclose(g, np.arange(6).reshape(2, 3)
                               * 27.211386245988 / 0.529177210903,
                               rtol=1e-6)


# ----------------------------------------------------------------------
# finite-difference Hessian
# ----------------------------------------------------------------------
def spring_factory(k, r0):
    def factory(directory, displacement):
        return Springs(k, r0)
    return factory


@pytest.mark.parametrize("scheme,tol", [("central", 1e-4), ("forward", 2e-2)])
def test_fd_hessian_matches_analytic(tmp_path, scheme, tol):
    atoms, k, r0 = water_like()
    # a slightly strained geometry so that forces are non-zero
    atoms.positions[1] += [0.02, 0.01, 0.0]
    fd = FiniteDifferenceHessian(atoms, spring_factory(k, r0),
                                 tmp_path / "h", step=0.002, scheme=scheme)
    assert fd.run() == len(fd.displacements())
    h_ref = spring_hessian(atoms.get_positions(), k, r0)
    h = fd.hessian()
    assert np.abs(h - h_ref).max() < tol * np.abs(h_ref).max()


def test_fd_hessian_resume_parts_and_metadata(tmp_path):
    atoms, k, r0 = water_like()
    fd = FiniteDifferenceHessian(atoms, spring_factory(k, r0), tmp_path,
                                 step=0.003)
    njobs = len(fd.displacements())
    assert njobs == 1 + 6 * len(atoms)
    n0 = fd.run(part=0, nparts=3)
    n1 = fd.run(part=1, nparts=3)
    assert len(fd.status()["pending"]) == njobs - n0 - n1
    # reopen and finish; nothing is recomputed
    fd2 = FiniteDifferenceHessian.open(tmp_path, spring_factory(k, r0))
    fd2.run(part=2, nparts=3)
    assert fd2.run() == 0
    assert len(fd2.status()["done"]) == njobs
    # a stale lock is skipped and can be cleared
    (tmp_path / "a0000x+" / "result.npz").unlink()
    (tmp_path / "a0000x+" / ".lock").write_text("dead")
    assert fd2.run() == 0
    assert fd2.clear_locks() == 1
    assert fd2.run() == 1
    # same workdir, different step: refused
    with pytest.raises(ValueError):
        FiniteDifferenceHessian(atoms, spring_factory(k, r0), tmp_path,
                                step=0.004)


def test_partial_hessian(tmp_path):
    atoms, k, r0 = water_like()
    fd = FiniteDifferenceHessian(atoms, spring_factory(k, r0), tmp_path,
                                 indices=[0, 1, 2])
    fd.run()
    h_ref = spring_hessian(atoms.get_positions(), k, r0)[:9, :9]
    np.testing.assert_allclose(fd.hessian(), h_ref, atol=1e-3)
    assert fd.normal_modes().nmodes == 3


# ----------------------------------------------------------------------
# line shapes
# ----------------------------------------------------------------------
def one_mode(w, s):
    return VibronicCouplings(np.array([w]), np.array([s * w]),
                             np.array([0]))


def peak_weights(ls, e00, w, kmax):
    out = []
    for kk in range(kmax):
        sel = np.abs(ls.energy - (e00 + kk * w)) < 0.5 * w
        out.append(ls.intensity[sel].sum() * ls.step)
    return np.array(out)


def test_zero_temperature_poisson_progression():
    w, s = 0.18, 1.3
    ls = vibronic_lineshape(one_mode(w, s), 0.0, "absorption", sigma=0.01,
                            step=5e-4)
    e00 = -s * w
    expected = [math.exp(-s) * s**kk / math.factorial(kk) for kk in range(6)]
    np.testing.assert_allclose(peak_weights(ls, e00, w, 6), expected,
                               atol=2e-4)
    assert abs(ls.zero_zero - e00) < 1e-12


def test_normalization_mean_and_variance():
    rng = np.random.default_rng(2)
    w = rng.uniform(0.01, 0.2, 12)
    lam = rng.uniform(0.001, 0.03, 12)
    c = VibronicCouplings(w, lam, np.arange(12))
    t, lcl, sig = 300.0, 0.02, 0.01
    ls = vibronic_lineshape(c, t, "absorption", lambda_classical=lcl,
                            sigma=sig)
    kt = KB_EV * t
    var = (lam * w / np.tanh(w / (2 * kt))).sum() + 2 * lcl * kt + sig**2
    norm = ls.intensity.sum() * ls.step
    mean, std = ls.moments()
    assert abs(norm - 1) < 1e-10
    assert abs(mean) < 1e-8
    assert abs(std**2 - var) < 1e-6 * var


def test_high_temperature_limit_is_marcus():
    # a soft mode at 2000 K behaves classically
    w, lam, t = 0.002, 0.15, 2000.0
    ls_q = vibronic_lineshape(one_mode(w, lam / w), t, "absorption",
                              step=2e-3)
    ls_c = vibronic_lineshape(None, t, "absorption", lambda_classical=lam,
                              step=2e-3)
    assert np.abs(ls_q.intensity - ls_c.intensity).max() < \
        2e-3 * ls_c.intensity.max()
    # classical_below moves the mode to the classical part
    ls_m = vibronic_lineshape(one_mode(w, lam / w), t, "absorption",
                              classical_below=0.01,
                              step=2e-3)
    np.testing.assert_allclose(ls_m.intensity, ls_c.intensity, atol=1e-10)


def test_detailed_balance_and_processes():
    # no Gaussian/Lorentzian broadening (they break detailed balance); the
    # classical part converges the time integral and obeys it
    rng = np.random.default_rng(5)
    c = VibronicCouplings(rng.uniform(0.05, 0.2, 5),
                          rng.uniform(0.01, 0.05, 5), np.arange(5))
    t = 300.0
    kt = KB_EV * t
    ls = vibronic_lineshape(c, t, "absorption", lambda_classical=0.05,
                            step=1e-3)
    e00 = ls.zero_zero
    x = np.linspace(0.02, 0.25, 12)
    lo, hi = ls(e00 - x), ls(e00 + x)
    np.testing.assert_allclose(lo, hi * np.exp(-x / kt), rtol=2e-3)
    y = np.linspace(-0.4, 0.4, 41)
    for process, above in PROCESSES.items():
        other = vibronic_lineshape(c, t, process, lambda_classical=0.05,
                                   step=1e-3)
        ref = ls(y) if above else ls(-y)
        np.testing.assert_allclose(other(y), ref,
                                   atol=1e-9 * ls.intensity.max())
        assert abs(other.zero_zero - (e00 if above else -e00)) < 1e-12
    with pytest.raises(ValueError, match="unknown process"):
        vibronic_lineshape(c, t, "mirror", lambda_classical=0.05)


def test_unconverged_lineshape_warns():
    with pytest.warns(UserWarning, match="not converged"):
        vibronic_lineshape(one_mode(0.1, 1.0), 0.0, "absorption", gamma=1e-5,
                           step=1e-3)


def test_ensemble_spectrum():
    c = one_mode(0.15, 0.8)
    ls = vibronic_lineshape(c, 300.0, "absorption", lambda_classical=0.02,
                            sigma=0.02)
    sites = np.array([-1.2, -1.0, -0.7])
    w = np.array([1.0, 0.5, 2.0])
    e, inten = ensemble_spectrum(sites, ls, weights=w)
    assert abs(np.trapezoid(inten, e) - w.sum()) < 1e-3
    # single site: the line shape itself, shifted (sites between grid
    # points are split linearly: error O(step^2))
    e1 = np.linspace(-1.6, -0.6, 501)
    _, i1 = ensemble_spectrum([-1.0], ls, energies=e1)
    np.testing.assert_allclose(i1, ls(e1 + 1.0), atol=1e-3 * i1.max())
    # off-grid site energy: linear weights keep the first moment
    e3 = np.linspace(-2.5, 0.5, 4001)
    _, i1 = ensemble_spectrum([-1.0], ls, energies=e3)
    _, i2 = ensemble_spectrum([-1.0 + 0.3 * ls.step], ls, energies=e3)
    m1, m2 = (np.trapezoid(i * e3, e3) / np.trapezoid(i, e3)
              for i in (i1, i2))
    assert abs(m1 + 1.0) < 1e-5
    assert abs(m2 - m1 - 0.3 * ls.step) < 1e-6


def test_gap_gradient_from_orb(tmp_path):
    import h5py
    f_gs = np.array([[0.01, 0.0, -0.02], [0.0, 0.03, 0.0]])
    f_es = np.array([[0.05, -0.01, 0.0], [0.02, 0.0, 0.01]])
    with h5py.File(tmp_path / "x.orb", "w") as f:
        q = f.create_group("QMdata")
        q["forces"] = f_gs
        q["excited_state_forces"] = f_es
    g = read_gap_gradient(tmp_path / "x.orb")
    np.testing.assert_allclose(g, (f_gs - f_es) * 27.211386245988
                               / 0.529177210903, rtol=1e-6)
    np.save(tmp_path / "g.npy", g)
    np.testing.assert_allclose(read_gap_gradient(tmp_path / "g.npy"), g)


@pytest.mark.skipif(not os.environ.get("PYXTP_RUN_DFTGWBSE"),
                    reason="runs xtp DFT; set PYXTP_RUN_DFTGWBSE=1")
def test_xtp_water_hessian(tmp_path):
    from ase.build import molecule
    from ase.optimize import BFGS
    from pyxtp.vibronic import xtp_ground_state_factory

    settings = {"dftpackage/basisset": "def2-svp",
                "dftpackage/auxbasisset": "aux-def2-svp",
                "dftpackage/functional": "XC_HYB_GGA_XC_PBEH"}
    factory = xtp_ground_state_factory(settings)
    atoms = molecule("H2O")
    optdir = tmp_path / "opt"
    optdir.mkdir()
    cwd = os.getcwd()
    os.chdir(optdir)
    try:
        atoms.calc = factory(optdir, None)
        BFGS(atoms, logfile=None).run(fmax=1e-3)
    finally:
        os.chdir(cwd)
    atoms.calc = None
    fd = FiniteDifferenceHessian(atoms, factory, tmp_path / "hess",
                                 step=0.005)
    fd.run()
    modes = fd.normal_modes()
    assert np.abs(fd.reference_forces()).max() < 2e-3
    assert fd.asymmetry < 0.05
    wn = modes.wavenumbers
    assert modes.nmodes == 3 and (wn > 0).all()
    assert 1550 < wn[0] < 1700 and 3700 < wn[1] < 4000 and 3800 < wn[2] < 4100
    # the displaced SCFs were warm-started from the reference orbitals
    log = next((tmp_path / "hess" / "a0000x+").glob("*.out")).read_text()
    assert "orbfile" in log.lower() or "guess" in log.lower()
