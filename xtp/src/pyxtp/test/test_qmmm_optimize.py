"""QM/MM geometry optimizer (pyxtp.qmmm_optimize) against a fake
xtp_parallel: a harmonic model of the QM region that reads the options and
job file like the real program and writes its job output in the same
format. Tests the plumbing (step directories, qmregion geometry/guess,
units, parsing), not the physics."""
import os
import stat
import sys
import textwrap
import unittest

import numpy as np
import pytest

try:
    from ase import units
    from pyxtp.qmmm_optimize import QMMMCalculator, QMMMJobError, optimize
except ImportError:
    raise unittest.SkipTest('pyxtp not available')

# the mapped geometry (Angstrom) and the minimum of the model, rigidly moved
# away from it: the optimizer has to translate the molecule, which the
# segment mapping would undo
MAPPED = [[0.0, 0.0, 0.0], [1.09, 0.0, 0.0], [-0.36, 1.03, 0.0]]
MINIMUM = [[0.30, -0.20, 0.10], [1.25, -0.20, 0.10], [0.00, 0.80, 0.10]]
# environment atoms of the fake checkpoint (Angstrom): a water and a
# virtual site
ENVIRONMENT = [[3.2, 0.4, 0.3], [3.8, 1.1, 0.3], [3.5, -0.4, 0.8],
               [3.3, 0.3, 0.4]]

FAKE = textwrap.dedent('''\
    #!{python}
    # fake xtp_parallel -e qmmm: E = k/2 |r - r0|^2 (Hartree, bohr)
    import os, sys, xml.etree.ElementTree as ET
    import numpy as np
    bohr = 0.529177210903
    k = 0.5
    els = ["C", "H", "H"]
    r0 = np.array({minimum}) / bohr
    opts = ET.parse("qmmm.xml").getroot().find("qmmm/regions/qmregion")
    assert opts.findtext("compute_forces") == "true"
    geo = opts.findtext("geometry")
    if geo:
        rows = [l.split() for l in open(geo).read().splitlines()[2:] if l.strip()]
        assert [r[0] for r in rows] == els
        r = np.array([[float(v) for v in row[1:]] for row in rows]) / bohr
    else:
        r = np.array({mapped}) / bohr
    guess = opts.findtext("guess")
    if guess:
        assert os.path.isfile(guess), guess
    with open("guess_used.txt", "w") as f:
        f.write(guess or "")
    polar = [p.findtext("guess") or "" for p in
             ET.parse("qmmm.xml").getroot().findall("qmmm/regions/polarregion")]
    with open("polar_guess_used.txt", "w") as f:
        f.write("\\n".join(polar))
    e = 0.5 * k * ((r - r0) ** 2).sum() - 40.0
    f = -k * (r - r0)
    tree = ET.parse("jobs.xml")
    job = tree.getroot().find("job")
    job.find("status").text = "COMPLETE"
    out = ET.SubElement(job, "output")
    regs = ET.SubElement(out, "regions")
    reg = ET.SubElement(regs, "region", type="qmregion", id="0")
    ET.SubElement(reg, "Initial_State").text = "s1"
    ET.SubElement(reg, "Final_State").text = "s1"
    fo = ET.SubElement(reg, "forces", state="s1", unit="Hartree/Bohr")
    go = ET.SubElement(reg, "geometry", unit="Angstrom")
    for i, el in enumerate(els):
        ET.SubElement(fo, "atom", element=el, id=str(i)).text = " ".join(
            "%+1.10e" % v for v in f[i])
        ET.SubElement(go, "atom", element=el, id=str(i)).text = " ".join(
            "%+1.10f" % v for v in r[i] * bohr)
    ET.SubElement(out, "E_tot").text = "%.10f" % (e * 27.211386245988)
    tree.write("jobs.xml")
    d = os.path.join("QMMM", "frame_0", "job_0")
    os.makedirs(d)
    open(os.path.join(d, "qmregion_final.orb"), "w").write("orbitals")
    for n in (2, 10, 9):
        open(os.path.join(d, "checkpoint_iter_%d.hdf5" % n), "w").write("cpt")
    if os.environ.get("FAKE_H5"):
        # a real checkpoint: QM atoms and a polar region with three atoms and
        # one virtual site, positions in bohr as qmmm writes them
        import h5py
        sdt = np.dtype([("index", "<i8"), ("type", h5py.string_dtype()),
                        ("posX", "<f8"), ("posY", "<f8"), ("posZ", "<f8")])
        with h5py.File(os.path.join(d, "checkpoint_iter_10.hdf5"), "w") as h:
            q = h.create_group("region_0")
            q.attrs["type"] = np.array(["qmregion"], dtype=h5py.string_dtype())
            qa = np.zeros((len(els), 1), dtype=[("index", "<i8"),
                          ("element", h5py.string_dtype()), ("posX", "<f8"),
                          ("posY", "<f8"), ("posZ", "<f8")])
            for i, el in enumerate(els):
                qa[i, 0] = (i, el, *r[i])
            q.create_dataset("QMdata/qmmolecule/qmatoms", data=qa)
            p = h.create_group("region_1")
            p.attrs["type"] = np.array(["polarregion"], dtype=h5py.string_dtype())
            sites = np.zeros((4, 1), dtype=sdt)
            env = np.array({environment}) / bohr
            for i, el in enumerate(["O", "H", "H", "X"]):
                sites[i, 0] = (i, el, *env[i])
            p.create_dataset("segments/PolarSegment_7/polarsites", data=sites)
    ''')


@pytest.fixture
def setup(tmp_path):
    fake = tmp_path / "xtp_parallel"
    fake.write_text(FAKE.format(python=sys.executable, minimum=MINIMUM,
                                mapped=MAPPED, environment=ENVIRONMENT))
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    s = tmp_path / "setup"
    s.mkdir()
    (s / "state.hdf5").write_text("state")
    (s / "qmmm.xml").write_text(
        "<options><qmmm><job_file>jobs.xml</job_file><regions><qmregion>"
        "<id>0</id><compute_forces>false</compute_forces></qmregion>"
        "<polarregion><id>1</id></polarregion>"
        "<polarregion><id>2</id><guess>stale.hdf5</guess></polarregion>"
        "</regions></qmmm></options>")
    (s / "jobs.xml").write_text(
        "<jobs><job><id>0</id><tag>mol_0:s1</tag><input/>"
        "<status>COMPLETE</status><output>old</output></job></jobs>")
    return tmp_path, str(fake)


def test_sella_keeps_rigid_body_motion(setup):
    # the model minimum is a rigid translation of the mapped geometry: an
    # optimizer that projects out translations (Sella's default for free
    # molecules) would stop at the mapped position with a net force
    pytest.importorskip("sella")
    tmp_path, fake = setup
    calc = QMMMCalculator(setup=str(tmp_path / "setup"), options="qmmm.xml",
                          workdir=str(tmp_path / "opt"), executable=fake)
    atoms, converged = optimize(calc, fmax=1e-3, steps=80, optimizer="Sella",
                                logfile=None)
    assert converged
    assert np.abs(atoms.get_positions() - np.array(MINIMUM)).max() < 1e-3
    assert np.abs(atoms.get_forces()).max() < 1e-3


def test_sella_without_rigid_body_motion(setup):
    # with rigid_body=False the centroid of the QM region stays where the
    # mapping put it, although the model pulls it away (net force left)
    pytest.importorskip("sella")
    tmp_path, fake = setup
    calc = QMMMCalculator(setup=str(tmp_path / "setup"), options="qmmm.xml",
                          workdir=str(tmp_path / "opt"), executable=fake)
    atoms, converged = optimize(calc, fmax=1e-3, steps=20, optimizer="Sella",
                                logfile=None, rigid_body=False)
    assert converged
    assert np.allclose(atoms.get_positions().mean(axis=0),
                       np.array(MAPPED).mean(axis=0), atol=1e-6)
    assert np.linalg.norm(atoms.get_forces().sum(axis=0)) > 1e-2
    with pytest.raises(ValueError, match="Sella"):
        optimize(calc, optimizer="BFGS", rigid_body=False, logfile=None)


def test_optimizes_to_the_model_minimum(setup):
    tmp_path, fake = setup
    calc = QMMMCalculator(setup=str(tmp_path / "setup"), options="qmmm.xml",
                          workdir=str(tmp_path / "opt"), executable=fake)
    atoms = calc.start_atoms()
    assert np.allclose(atoms.get_positions(), MAPPED)
    atoms, converged = optimize(calc, fmax=1e-3, steps=50, atoms=atoms,
                                logfile=None)
    assert converged
    # rigid translation kept: nothing re-aligns the molecule
    assert np.abs(atoms.get_positions() - np.array(MINIMUM)).max() < 1e-3
    # the first job's result is reused for the start geometry
    assert calc.history[0]["step"] == 0 and len(calc.history) == calc.nsteps
    # every step after the first runs at the given geometry with the
    # previous step's orbitals as guess
    for i in range(1, calc.nsteps):
        stepdir = tmp_path / "opt" / f"step_{i:03d}"
        assert (stepdir / "geometry.xyz").is_file()
        guess = (stepdir / "guess_used.txt").read_text()
        assert guess.endswith(os.path.join(f"step_{i - 1:03d}", "QMMM",
                                           "frame_0", "job_0",
                                           "qmregion_final.orb"))
        # every polar region starts from the last checkpoint of the
        # previous step (largest iteration number, not the last listed)
        polar = (stepdir / "polar_guess_used.txt").read_text().split("\n")
        expected = os.path.join(f"step_{i - 1:03d}", "QMMM", "frame_0",
                                "job_0", "checkpoint_iter_10.hdf5")
        assert len(polar) == 2 and all(p.endswith(expected) for p in polar)
    assert (tmp_path / "opt" / "step_000" / "guess_used.txt").read_text() == ""
    # no stale guess from the user's options in the first step
    assert (tmp_path / "opt" / "step_000" /
            "polar_guess_used.txt").read_text() == "\n"


def test_units(setup):
    tmp_path, fake = setup
    calc = QMMMCalculator(setup=str(tmp_path / "setup"), options="qmmm.xml",
                          workdir=str(tmp_path / "opt"), executable=fake)
    atoms = calc.start_atoms()
    f = atoms.get_forces()
    expected = -0.5 * (np.array(MAPPED) - np.array(MINIMUM)) / units.Bohr
    assert np.allclose(f, expected * units.Hartree / units.Bohr, rtol=1e-8)
    # finite difference of the energy along one coordinate
    h = 1e-4
    p = atoms.get_positions()
    e = []
    for s in (1, -1):
        q = p.copy()
        q[1, 0] += s * h
        atoms.set_positions(q)
        e.append(atoms.get_potential_energy())
    assert (e[0] - e[1]) / (2 * h) == pytest.approx(-f[1, 0], rel=1e-5)


def test_failed_job_is_reported(setup):
    tmp_path, _ = setup
    bad = tmp_path / "bad"
    bad.write_text(f"#!{sys.executable}\nimport sys\nsys.exit(3)\n")
    bad.chmod(bad.stat().st_mode | stat.S_IEXEC)
    calc = QMMMCalculator(setup=str(tmp_path / "setup"), options="qmmm.xml",
                          workdir=str(tmp_path / "opt"), executable=str(bad))
    with pytest.raises(QMMMJobError, match="exited with 3"):
        calc.start_atoms()


# ---- D3 dispersion ---------------------------------------------------------
def _d3():
    try:
        import dftd3.interface  # noqa: F401
    except ImportError:
        pytest.skip("dftd3 (simple-dftd3) not available")
    from pyxtp.qmmm_optimize import D3Dispersion
    return D3Dispersion


QM = (["C", "H", "H", "H", "O"],
      np.array([[0.0, 0.0, 0.0], [1.09, 0.0, 0.0], [-0.36, 1.03, 0.0],
                [-0.36, -0.51, 0.89], [-0.6, -0.4, -1.2]]))
ENV = (["O", "H", "H", "C", "H"],
       np.array([[3.0, 0.2, 0.1], [3.6, 0.9, 0.1], [3.3, -0.6, 0.6],
                 [-1.0, 3.5, 0.5], [-1.2, 4.5, 0.6]]))


def test_d3_qm_only_matches_dftd3_ase():
    D3Dispersion = _d3()
    from ase import Atoms
    from dftd3.ase import DFTD3
    d = D3Dispersion(method="pbe0", environment=False)
    e, f = d(*QM)
    atoms = Atoms(QM[0], positions=QM[1])
    atoms.calc = DFTD3(method="pbe0", damping="d3bj")
    assert e == pytest.approx(atoms.get_potential_energy(), abs=1e-10)
    assert np.allclose(f, atoms.get_forces(), atol=1e-10)


@pytest.mark.parametrize("atm", [False, True])
def test_d3_with_environment(atm):
    D3Dispersion = _d3()
    d = D3Dispersion(method="pbe0", atm=atm)
    d.set_environment(ENV[0], ENV[1], QM[1])
    e, f = d(*QM)
    # interaction + QM part: E(QM + env) - E(env)
    full = D3Dispersion(method="pbe0", atm=atm, environment=False)
    e_all = full(QM[0] + ENV[0], np.vstack([QM[1], ENV[1]]))[0]
    e_env = full(*ENV)[0]
    assert e == pytest.approx(e_all - e_env, abs=1e-10)
    assert e != pytest.approx(full(*QM)[0], abs=1e-6)
    # forces on the QM atoms with the environment frozen
    h = 1e-4
    for i, k in [(0, 0), (4, 2), (2, 1)]:
        p = QM[1].copy()
        p[i, k] += h
        ep = d(QM[0], p)[0]
        p[i, k] -= 2 * h
        em = d(QM[0], p)[0]
        assert -(ep - em) / (2 * h) == pytest.approx(f[i, k], abs=1e-8)


def test_d3_environment_cutoff():
    D3Dispersion = _d3()
    d = D3Dispersion(method="pbe0", cutoff=2.6)
    d.set_environment(ENV[0], ENV[1], QM[1])
    assert d.env_symbols == ["O", "H"]  # the water's far H and the CH out


def test_optimizer_with_d3(setup, monkeypatch):
    D3Dispersion = _d3()
    from pyxtp.qmmm_optimize import environment_from_checkpoint
    tmp_path, fake = setup
    monkeypatch.setenv("FAKE_H5", "1")
    disp = D3Dispersion(method="pbe0")
    calc = QMMMCalculator(setup=str(tmp_path / "setup"), options="qmmm.xml",
                          workdir=str(tmp_path / "opt"), executable=fake,
                          terms=[disp])
    atoms = calc.start_atoms()
    # environment from the first checkpoint, without the virtual site
    assert disp.env_symbols == ["O", "H", "H"]
    assert np.allclose(disp.env_positions, np.array(ENVIRONMENT[:3]))
    cpt = os.path.join(calc.history[0]["stepdir"], "QMMM", "frame_0", "job_0",
                       "checkpoint_iter_10.hdf5")
    env_s, _, qm_s, qm_p = environment_from_checkpoint(cpt)
    assert env_s == ["O", "H", "H"] and qm_s == ["C", "H", "H"]
    assert np.allclose(qm_p, MAPPED)
    # energy and forces are those of the model plus D3
    r0 = calc.history[0]
    e_d3, f_d3 = disp(r0["elements"], r0["positions"])
    assert r0["energy_d3"] == pytest.approx(e_d3)
    assert r0["energy"] == pytest.approx(r0["energy_qmmm"] + e_d3)
    model = -0.5 * (np.array(MAPPED) - np.array(MINIMUM)) / units.Bohr
    assert np.allclose(atoms.get_forces(),
                       model * units.Hartree / units.Bohr + f_d3, atol=1e-8)
    atoms, converged = optimize(calc, fmax=1e-3, steps=60, atoms=atoms,
                                logfile=None)
    assert converged
    # D3 shifts the minimum of the model slightly
    shift = np.abs(atoms.get_positions() - np.array(MINIMUM)).max()
    assert 1e-4 < shift < 0.05


def test_lennard_jones():
    from pyxtp.qmmm_optimize import LennardJones
    qm = {"C": (3.4, 0.36), "H": (2.6, 0.06), "O": (3.07, 0.88)}
    env = {"O": (3.15, 0.64), "C": (3.4, 0.36)}  # environment H: no LJ
    lj = LennardJones(qm, env, cutoff=None)
    lj.set_environment(ENV[0], ENV[1], QM[1])
    assert lj.env_symbols == ["O", "C"] and lj.missing == ["H"]
    e, f = lj(*QM)
    # direct sum, Lorentz-Berthelot, kJ/mol -> eV
    ref = 0.0
    for s, p in zip(QM[0], QM[1]):
        for t, q in zip(ENV[0], ENV[1]):
            if t not in env:
                continue
            sig = 0.5 * (qm[s][0] + env[t][0])
            eps = (qm[s][1] * env[t][1]) ** 0.5
            r = np.linalg.norm(p - q)
            ref += 4 * eps * ((sig / r) ** 12 - (sig / r) ** 6)
    assert e == pytest.approx(ref * units.kJ / units.mol, rel=1e-12)
    h = 1e-5
    for i, k in [(0, 0), (4, 1), (1, 2)]:
        p = QM[1].copy()
        p[i, k] += h
        ep = lj(QM[0], p)[0]
        p[i, k] -= 2 * h
        em = lj(QM[0], p)[0]
        assert -(ep - em) / (2 * h) == pytest.approx(f[i, k], rel=1e-6, abs=1e-10)
    # per-atom parameters give the same as the element table
    lj2 = LennardJones([qm[s] for s in QM[0]], env, cutoff=None)
    lj2.set_environment(ENV[0], ENV[1], QM[1])
    assert lj2(*QM)[0] == pytest.approx(e, rel=1e-14)


def test_optimizer_with_lj_and_d3(setup, monkeypatch):
    D3Dispersion = _d3()
    from pyxtp.qmmm_optimize import LennardJones
    tmp_path, fake = setup
    monkeypatch.setenv("FAKE_H5", "1")
    lj = LennardJones({"C": (3.4, 0.36), "H": (2.6, 0.06)},
                      {"O": (3.15, 0.64), "H": (0.0, 0.0)})
    d3 = D3Dispersion(method="pbe0", environment=False)
    calc = QMMMCalculator(setup=str(tmp_path / "setup"), options="qmmm.xml",
                          workdir=str(tmp_path / "opt"), executable=fake,
                          terms=[lj, d3])
    atoms = calc.start_atoms()
    r0 = calc.history[0]
    assert r0["energy"] == pytest.approx(
        r0["energy_qmmm"] + r0["energy_lj"] + r0["energy_d3"])
    atoms, converged = optimize(calc, fmax=1e-3, steps=60, atoms=atoms,
                                logfile=None)
    assert converged
    with pytest.raises(ValueError, match="same name"):
        QMMMCalculator(setup=str(tmp_path / "setup"), options="qmmm.xml",
                       workdir=str(tmp_path / "opt2"), executable=fake,
                       terms=[d3, d3])
