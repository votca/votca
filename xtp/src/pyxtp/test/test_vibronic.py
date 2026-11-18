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
    from pyxtp.vibronic import QPLevels, Spectrum, level_spectrum
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


# ----------------------------------------------------------------------
# spectra of several levels
# ----------------------------------------------------------------------
def toy_levels():
    """Spring molecule with four 'QP levels' (HOMO-1 .. LUMO+1)."""
    atoms, k, r0 = water_like()
    modes = NormalModes.from_hessian(
        atoms, spring_hessian(atoms.get_positions(), k, r0))
    rng = np.random.default_rng(11)
    g = rng.normal(scale=0.4, size=(4, len(atoms), 3))
    levels = QPLevels(levels=np.array([8, 9, 10, 11]), homo=9,
                      energies=np.array([-8.0, -6.0, -1.0, 0.5]),
                      gradients=g)
    return modes, levels


def test_level_spectrum_selects_levels_and_axis():
    modes, levels = toy_levels()
    kw = dict(temperature=300.0, sigma=0.03)
    add = level_spectrum(modes, levels, "electron-addition", **kw)
    assert list(add.contributions) == ["LUMO", "LUMO+1"]
    rem = level_spectrum(modes, levels, "electron-removal", **kw)
    assert list(rem.contributions) == ["HOMO-1", "HOMO"]
    ion = level_spectrum(modes, levels, "ionization", **kw)
    # every level integrates to one; the mean is the vertical energy
    for sp, label, evert in ((add, "LUMO", -1.0), (rem, "HOMO", -6.0),
                             (ion, "HOMO", 6.0), (ion, "HOMO-1", 8.0)):
        y = sp.contributions[label]
        assert abs(np.trapezoid(y, sp.energy) - 1.0) < 2e-3
        assert abs(np.trapezoid(y * sp.energy, sp.energy) - evert) < 5e-3
    # the IP axis is the mirror image of the level axis
    np.testing.assert_allclose(
        np.interp(-rem.energy, ion.energy, ion.contributions["HOMO"]),
        rem.contributions["HOMO"], atol=2e-3 * rem.total.max())
    only = level_spectrum(modes, levels, "electron-addition",
                          select=["lumo+1"], **kw)
    assert list(only.contributions) == ["LUMO+1"]
    with pytest.raises(ValueError):
        level_spectrum(modes, levels, "absorption", **kw)


def test_level_spectrum_sites_and_rigid_shift():
    modes, levels = toy_levels()
    rng = np.random.default_rng(4)
    lumo_sites = -1.3 + 0.1 * rng.standard_normal(200)
    sp = level_spectrum(modes, levels, "electron-addition", 300.0,
                        site_energies={"LUMO": lumo_sites}, sigma=0.03)
    assert "200 site energies" in sp.notes["LUMO"]
    assert "shifted by +1.5000 eV" in sp.notes["LUMO+1"]
    for label, mean in (("LUMO", lumo_sites.mean()),
                        ("LUMO+1", lumo_sites.mean() + 1.5)):
        y = sp.contributions[label]
        area = np.trapezoid(y, sp.energy)
        assert abs(area - 200) < 0.5
        assert abs(np.trapezoid(y * sp.energy, sp.energy) / area - mean) \
            < 5e-3
    with pytest.raises(ValueError, match="without gradients"):
        level_spectrum(modes, levels, "electron-addition", 300.0,
                       site_energies={"homo": lumo_sites}, sigma=0.03)


def test_spectrum_broadening_onset_and_file(tmp_path):
    e = np.linspace(-3.0, 3.0, 3001)
    s0, mu = 0.2, 0.4
    gauss = np.exp(-0.5 * ((e - mu) / s0) ** 2)
    sp = Spectrum(e, {"LUMO": gauss, "LUMO+1": 0.0 * gauss},
                  "electron-addition", {"LUMO": 0.1, "LUMO+1": 0.2})
    # tangent at the inflection point of a Gaussian hits zero at mu - 2 s
    res = sp.onset()
    assert abs(res["onset"] - (mu - 2 * s0)) < 2e-3
    assert abs(res["peak"] - mu) < 2e-3
    high = Spectrum(e, {"HOMO": gauss}, "electron-removal").onset()
    assert abs(high["onset"] - (mu + 2 * s0)) < 2e-3
    # Gaussian broadening adds the variances and keeps the area
    b = sp.broadened(0.15)
    area = np.trapezoid(b.total, b.energy)
    assert abs(area - np.trapezoid(gauss, e)) < 1e-6
    var = np.trapezoid(b.total * (b.energy - mu) ** 2, b.energy) / area
    assert abs(var - (s0**2 + 0.15**2)) < 1e-4
    sp.save(tmp_path / "s.dat", header="test")
    back = Spectrum.load(tmp_path / "s.dat")
    assert list(back.contributions) == ["LUMO", "LUMO+1"]
    assert back.process == "electron-addition"
    assert back.reorganization == {"LUMO": 0.1, "LUMO+1": 0.2}
    np.testing.assert_allclose(back.total, sp.total)


def test_qp_levels_from_orb(tmp_path):
    import h5py
    g = np.arange(2 * 3 * 3, dtype=float).reshape(6, 3) * 1e-3
    with h5py.File(tmp_path / "x.orb", "w") as f:
        q = f.create_group("QMdata")
        q["qp_gradients"] = g
        q["qp_gradient_levels"] = np.array([4, 5])
        q["QPpert_energies"] = np.array([[-0.5], [-0.4], [0.1], [0.2]])
        q.attrs["occupied_levels"] = np.array([5])
        q.attrs["qpmin"] = np.array([3])
    lv = QPLevels.from_orb(tmp_path / "x.orb")
    assert lv.labels == ["HOMO", "LUMO"]
    hartree, bohr = 27.211386245988, 0.529177210903
    np.testing.assert_allclose(lv.energies, [-0.4 * hartree, 0.1 * hartree],
                               rtol=1e-6)
    np.testing.assert_allclose(lv.gradients[1], g[3:] * hartree / bohr,
                               rtol=1e-6)
    assert list(lv.select(True).levels) == [4]


def test_plots(tmp_path):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from pyxtp.vibronic.plot import (paper_style, plot_couplings,
                                     plot_spectrum)
    modes, levels = toy_levels()
    sp = level_spectrum(modes, levels, "electron-addition", 300.0,
                        sigma=0.03)
    with paper_style():
        fig = plot_spectrum(sp, broadening=0.2, onset=True,
                            experiment=(sp.energy, sp.total))
        assert fig.onset["onset"] < fig.onset["peak"]
        fig.savefig(tmp_path / "s.pdf")
        fig = plot_couplings(levels.couplings(modes)[0])
        fig.savefig(tmp_path / "c.pdf")
    assert (tmp_path / "s.pdf").stat().st_size > 1000
    assert (tmp_path / "c.pdf").stat().st_size > 1000


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


@pytest.mark.skipif(not os.environ.get("PYXTP_RUN_DFTGWBSE"),
                    reason="runs dftgwbse; set PYXTP_RUN_DFTGWBSE=1")
def test_xtp_qp_gradients(tmp_path):
    """Analytic QP gradients (G0W0) against finite differences."""
    import h5py
    from pyxtp import xtp

    pos = np.array([[0.0, 0.0, 0.119], [0.03, 0.743, -0.437],
                    [0.0, -0.763, -0.477]])
    direction = np.random.default_rng(1).normal(size=(3, 3))
    direction /= np.linalg.norm(direction)

    def run(positions, tag, gradients):
        workdir = tmp_path / tag
        workdir.mkdir()
        cwd = os.getcwd()
        os.chdir(workdir)
        try:
            calc = xtp()
            calc.set(**{"dftpackage/basisset": "def2-svp",
                        "dftpackage/auxbasisset": "aux-def2-svp",
                        "dftpackage/convergence_tightness": "verytight",
                        "gwbse/gw/mode": "G0W0", "gwbse/ranges": "explicit",
                        "gwbse/qpmin": "0", "gwbse/qpmax": "9",
                        "gwbse/rpamax": "23", "gwbse/bsemin": "0",
                        "gwbse/bsemax": "9", "gwbse/bse/exctotal": "2"})
            if gradients:
                calc.options.qp_gradients = "homo-1:lumo"
            atoms = Atoms("OH2", positions=positions)
            atoms.calc = calc
            atoms.get_potential_energy()
            orb = str(workdir / "H2O.orb")
            with h5py.File(orb, "r") as f:
                qp = np.array(f["QMdata"]["QPpert_energies"][()]).ravel()
            return orb, qp * 27.211386245988
        finally:
            os.chdir(cwd)

    orb, _ = run(pos, "ref", True)
    levels = QPLevels.from_orb(orb)
    assert levels.labels == ["HOMO-1", "HOMO", "LUMO"]
    h = 1e-3
    _, plus = run(pos + h * direction, "plus", False)
    _, minus = run(pos - h * direction, "minus", False)
    for i, level in enumerate(levels.levels):
        fd = (plus[level] - minus[level]) / (2 * h)
        analytic = (levels.gradients[i] * direction).sum()
        assert abs(analytic - fd) < 1e-4 * np.linalg.norm(levels.gradients[i])
        assert np.abs(levels.gradients[i].sum(axis=0)).max() < 1e-6

    # more levels later, from the orb file alone (tasks=gradients): the
    # stored ones are kept, the new ones equal those of a full run
    from pyxtp.vibronic import add_gradients
    settings = {"dftpackage/basisset": "def2-svp",
                "dftpackage/auxbasisset": "aux-def2-svp",
                "dftpackage/convergence_tightness": "verytight",
                "gwbse/gw/mode": "G0W0", "gwbse/ranges": "explicit",
                "gwbse/qpmin": "0", "gwbse/qpmax": "9",
                "gwbse/rpamax": "23", "gwbse/bsemin": "0",
                "gwbse/bsemax": "9", "gwbse/bse/exctotal": "2"}
    mtime = os.path.getmtime(tmp_path / "ref" / "H2O.out")
    add_gradients(orb, qp_levels="homo-2:lumo+1", settings=settings)
    more = QPLevels.from_orb(orb)
    assert more.labels == ["HOMO-2", "HOMO-1", "HOMO", "LUMO", "LUMO+1"]
    np.testing.assert_array_equal(more.gradients[1:4], levels.gradients)
    assert os.path.getmtime(tmp_path / "ref" / "H2O.out") == mtime
    log = (tmp_path / "ref" / "H2O_gradients.out").read_text()
    assert "3 requested level(s) are already stored" in log
    assert "DFT calculation" not in log and "GW" in log
    for i in (0, 4):
        level = more.levels[i]
        fd = (plus[level] - minus[level]) / (2 * h)
        analytic = (more.gradients[i] * direction).sum()
        assert abs(analytic - fd) < 1e-4 * np.linalg.norm(more.gradients[i])
    # other options than those of the stored run are refused
    wrong = dict(settings, **{"dftpackage/basisset": "def2-tzvp"})
    with pytest.raises(RuntimeError, match="basis set"):
        add_gradients(orb, qp_levels="lumo+2", settings=wrong)
