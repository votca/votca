"""Test the analytic BSE excited-state forces of the ASE calculator.

xtp(excited_state=...): dftgwbse runs of water, only with
PYXTP_RUN_DFTGWBSE=1.
"""
import os
import unittest

import numpy as np
import pytest
from ase import Atoms

try:
    from pyxtp import xtp
    from pyxtp.xtp import numeric_force_xtp
except ImportError:
    raise unittest.SkipTest('pyxtp not available')

RUN = os.environ.get("PYXTP_RUN_DFTGWBSE") == "1"


def make_calc(mode='G0W0', sigma='ppm'):
    calc = xtp(nthreads=1, excited_state='S1')
    calc.set(**{'dftpackage/basisset': 'def2-svp',
                'dftpackage/auxbasisset': 'aux-def2-svp',
                'dftpackage/functional': 'XC_HYB_GGA_XC_PBEH',
                'gwbse/gw/mode': mode,
                'gwbse/bse/useTDA': 'true',
                'gwbse/bse/exctotal': '3'})
    if mode == 'evGW':
        calc.set(**{'gwbse/gw/sc_limit': '1e-9',
                    'gwbse/gw/screening_update': 'every'})
    if sigma != 'ppm':
        # QP window below the levels whose QP solutions jump between
        # branches near this geometry with the full-frequency self-energy
        calc.set(**{'gwbse/gw/sigma_integrator': sigma,
                    'gwbse/ranges': 'explicit', 'gwbse/qpmin': '0',
                    'gwbse/qpmax': '9', 'gwbse/rpamax': '23',
                    'gwbse/bsemin': '0', 'gwbse/bsemax': '23'})
    return calc


def water(shift=0.0):
    return Atoms('OH2', positions=[[0, 0, 0.1 + shift], [0.96, 0.05, 0],
                                   [-0.24, 0.93, -0.07]])


@pytest.fixture(scope='module')
def derivatives(tmp_path_factory):
    """Skip when libint2 lacks the derivative integrals the forces need.

    The DFT engine then refuses compute_forces, before the SCF.
    """
    cwd = os.getcwd()
    os.chdir(tmp_path_factory.mktemp('derivative_probe'))
    try:
        calc = xtp(nthreads=1)
        calc.set(**{'dftpackage/basisset': 'def2-svp',
                    'dftpackage/auxbasisset': 'aux-def2-svp',
                    'dftpackage/xtpdft/compute_forces': 'true'})
        atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
        atoms.calc = calc
        atoms.get_potential_energy()
    except RuntimeError as error:
        if 'without derivative integrals' in str(error):
            pytest.skip('libint2 was built without derivative integrals')
        raise
    finally:
        os.chdir(cwd)


def test_state_label():
    assert xtp._parse_state(' T2 ') == 't2'
    with pytest.raises(ValueError):
        xtp._parse_state('q1')


@pytest.mark.skipif(not RUN, reason="set PYXTP_RUN_DFTGWBSE=1")
def test_excited_state_energy_and_forces(derivatives):
    atoms = water()
    atoms.calc = make_calc()
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    calc = atoms.calc
    assert forces.shape == (3, 3)
    # total energy of S1 = DFT + Omega(S1)
    assert energy == pytest.approx(calc.DFTenergy + calc.bse_singlet_energies[0])
    # forces sum to zero (translational invariance)
    assert np.abs(forces.sum(axis=0)).max() < 1e-4
    # differ from the ground-state forces
    assert np.abs(forces - calc.dft_forces).max() > 1e-2

    # moved atoms: the forces belong to the new geometry, not the old one
    atoms.positions[0, 2] += 0.05
    moved = atoms.get_forces()
    fresh = water(0.05)
    fresh.calc = make_calc()
    np.testing.assert_allclose(moved, fresh.get_forces(), atol=1e-6)
    assert np.abs(moved - forces).max() > 1e-3


@pytest.mark.skipif(not RUN, reason="set PYXTP_RUN_DFTGWBSE=1")
@pytest.mark.parametrize("mode,sigma", [("G0W0", "ppm"), ("evGW", "ppm"),
                                        ("evGW", "exact")])
def test_excited_state_forces_finite_difference(mode, sigma, derivatives):
    """Compare the analytic forces with finite differences.

    The analytic forces include the response of the QP corrections, so
    they are the derivative of the total S1 energy with GW recomputed at
    every geometry.
    """
    atoms = water()
    atoms.calc = make_calc(mode, sigma)
    force = atoms.get_forces()[1, 0]
    h = 5e-4  # Angstrom
    energies = []
    for step in (h, -h):
        displaced = water()
        displaced.positions[1, 0] += step
        displaced.calc = make_calc(mode, sigma)
        energies.append(displaced.get_potential_energy())
    fd = -(energies[0] - energies[1]) / (2 * h)
    assert force == pytest.approx(fd, rel=5e-4)


@pytest.mark.skipif(not RUN, reason="set PYXTP_RUN_DFTGWBSE=1")
def test_root_following(derivatives):
    """Follow the state from geometry to geometry with follow_root.

    From the second geometry on, the state is the one with the largest
    exciton overlap with the state of the previous geometry; for small
    steps it keeps its label, and energy and forces belong to it.
    """
    calc = make_calc()
    calc.follow_root = True
    calc.excited_state = calc.tracked_state = 's2'
    atoms = water()
    atoms.calc = calc
    atoms.get_forces()
    assert calc.root_reference is not None
    assert os.path.exists(calc.root_reference)
    for shift in (0.01, 0.02):
        atoms.positions = water(shift).positions
        forces = atoms.get_forces()
        assert calc.tracked_state == 's2'
        assert calc.excited_forces_state == 's2'
        assert atoms.get_potential_energy() == pytest.approx(
            calc.DFTenergy + calc.bse_singlet_energies[1])
        assert np.abs(forces.sum(axis=0)).max() < 1e-4
    with open(atoms.get_chemical_formula() + '.out') as out:
        assert 'Root following: s2' in out.read()


@pytest.mark.skipif(not RUN, reason="set PYXTP_RUN_DFTGWBSE=1")
def test_calculator_attached_to_new_atoms(derivatives):
    """Attach the calculator to new atoms.

    A scan with new Atoms per point must not return the results of the
    previous positions.
    """
    calc = make_calc()
    first = water()
    first.calc = calc
    f1 = first.get_forces()
    second = water(0.05)
    second.calc = calc
    f2 = second.get_forces()
    assert np.abs(f1 - f2).max() > 1e-3


@pytest.mark.skipif(not RUN, reason="set PYXTP_RUN_DFTGWBSE=1")
def test_numerical_ground_state_forces(derivatives):
    """Check the finite-difference ground-state forces.

    Without dftpackage.xtpdft.compute_forces the ground-state forces are
    finite differences; they must see the displaced positions (they were
    zero when the calculator kept a copy of the atoms) and agree with the
    analytic DFT forces. DFT only and two force components, as
    calculate_numerical_forces computes them (the full set took most of
    the time of this test in the CI Debug build).
    """
    dft_only = {'tasks': 'input,dft,parse'}
    analytic = make_calc()
    analytic.excited_state = analytic.tracked_state = None
    analytic.set(**dft_only, **{'dftpackage/xtpdft/compute_forces': 'true'})
    atoms = water()
    atoms.calc = analytic
    reference = atoms.get_forces()
    assert np.abs(reference).max() > 1e-2
    numerical = make_calc()
    numerical.excited_state = numerical.tracked_state = None
    numerical.set(**dft_only)
    for a, i in ((1, 0), (0, 2)):
        other = water()
        other.calc = xtp(options=numerical.options, nthreads=1)
        force = numeric_force_xtp(other, a, i, 1e-3,
                                  numerical.option_forces)
        assert force == pytest.approx(reference[a, i], abs=2e-3)
