import os
import unittest

import numpy as np
import pytest

try:
    from pyxtp.finite_field import (FiniteFieldResult, OrbData,
                                    central_first, central_second,
                                    match_states, mixing_groups, mo_rotation,
                                    state_overlaps, write_xml)
except ImportError:
    raise unittest.SkipTest('pyxtp not available')


def orbdata(X, Y, mos, vmin=0, homo=1, cmax=4):
    n = X.shape[0]
    return OrbData(0.0, np.arange(n, dtype=float), X, Y, mos, vmin, homo,
                   cmax, ["H"], np.zeros((1, 3)), np.ones(1))


def random_states(rng, n=4, nv=2, nc=3, tda=False):
    X = rng.normal(size=(n, nv, nc))
    Y = None if tda else 0.1 * rng.normal(size=(n, nv, nc))
    return X, Y


def test_central_differences_exact_for_polynomials():
    h = 1e-2
    f = {k: 2.0 + 3.0 * (k * h) - 0.5 * 7.0 * (k * h) ** 2
         for k in (-2, -1, 0, 1, 2)}
    for order in (2, 4):
        assert central_first(f, h, order) == pytest.approx(3.0, rel=1e-10)
        assert central_second(f, h, order) == pytest.approx(-7.0, rel=1e-8)
    # the 4th-order formula is exact for a cubic term, the 2nd-order not
    g = {k: (k * h) ** 3 for k in (-2, -1, 0, 1, 2)}
    assert central_first(g, h, 4) == pytest.approx(0.0, abs=1e-12)
    assert central_first(g, h, 2) == pytest.approx(h ** 2, rel=1e-10)


def test_mo_rotation_recovers_rotation():
    rng = np.random.default_rng(1)
    C = rng.normal(size=(5, 5))
    U = np.linalg.qr(rng.normal(size=(5, 5)))[0]
    assert np.allclose(mo_rotation(C, C @ U), U)


def test_state_overlaps_with_permutation_and_mo_signs():
    rng = np.random.default_rng(2)
    X, Y = random_states(rng)
    mos = rng.normal(size=(5, 5))
    ref = orbdata(X, Y, mos)
    # same states, reordered, with the signs of MOs 0 and 3 flipped
    perm = [2, 0, 3, 1]
    signs = np.array([-1, 1, 1, -1, 1.0])
    occ, vir = signs[:2], signs[2:]
    flip = lambda A: A * occ[None, :, None] * vir[None, None, :]
    new = orbdata(flip(X[perm]), flip(Y[perm]), mos * signs[None, :])
    ov = state_overlaps(ref, new)
    cols, vals = match_states(ov, 4)
    assert list(cols) == [perm.index(n) for n in range(4)]
    assert np.allclose(vals, 1.0)


def test_state_overlaps_tda():
    rng = np.random.default_rng(3)
    X, _ = random_states(rng, tda=True)
    mos = rng.normal(size=(5, 5))
    ov = state_overlaps(orbdata(X, None, mos), orbdata(X, None, mos))
    assert np.allclose(np.diag(ov), 1.0)


def test_mixing_groups():
    # field states: 0 clean, 1 and 2 mixed 50/50
    ov = np.array([[1.0, 0.0, 0.0],
                   [0.0, 0.7, 0.7],
                   [0.0, 0.7, 0.7]])
    assert mixing_groups([np.eye(3), ov], 3, 0.9) == [[0], [1, 2]]
    assert mixing_groups([np.eye(3)], 3, 0.9) == [[0], [1], [2]]


def test_xml_output(tmp_path):
    # parses the file written just above, not untrusted input
    import xml.etree.ElementTree as ET  # nosec B405

    nan = np.nan
    r = FiniteFieldResult(
        "z", 1e-3, 2, "singlet", np.array([0.1, 0.2]),
        np.array([nan, nan, 0.9]), np.array([nan, nan, 20.0]),
        np.array([[nan, nan, -0.4], [nan, nan, 1.2]]),
        np.array([[nan, nan, 1.0], [nan, nan, 2.0]]),
        np.array([1.0, 0.5]), [[0], [1]], [np.array([nan, nan, -0.4]),
                                           np.array([nan, nan, 1.2])])
    write_xml(r, tmp_path / "ff.xml")
    root = ET.parse(tmp_path / "ff.xml").getroot()  # nosec B314
    assert root.tag == "finite_field" and root.get("components") == "z"
    assert float(root.find("ground_state/dipole").get("z")) == pytest.approx(0.9)
    assert root.find("ground_state/dipole").get("x") is None
    states = root.findall("state")
    assert [s.get("reliable") for s in states] == ["true", "false"]
    assert float(states[1].find("dipole").get("z")) == pytest.approx(2.1)


@pytest.mark.skipif(not os.environ.get("PYXTP_RUN_DFTGWBSE"),
                    reason="runs dftgwbse; set PYXTP_RUN_DFTGWBSE=1")
def test_carbon_monoxide(tmp_path):
    from ase.build import molecule
    from pyxtp import xtp
    from pyxtp.finite_field import FiniteFieldDipoles

    calc = xtp()
    calc.set(basis="def2-svp", **{"dftpackage/auxbasisset": "aux-def2-svp",
                                 "gwbse/gw/mode": "G0W0"})
    atoms = molecule("CO", positions=([0, 0, 0], [1.13, 0, 0]))
    ff = FiniteFieldDipoles(atoms, options=calc.options, nstates=2,
                            components="x", workdir=tmp_path)
    result = ff.run()
    assert np.isfinite(result.ground_dipole[0])
    assert np.all(np.isfinite(result.delta_dipole[:, 0]))
    assert np.all(result.min_overlap > 0.99)
    # only x was perturbed
    assert np.all(np.isnan(result.delta_dipole[:, 1:]))
