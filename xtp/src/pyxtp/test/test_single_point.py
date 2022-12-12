def test_single_point():
    import unittest

    from ase.build import molecule
    try:
        from pyxtp import xtp 
    except ImportError:
        raise unittest.SkipTest('xtp not available')

    
    # calculate H2
    calculator = xtp()
    atoms = molecule('H2', calculator=calculator)
    atoms.center(vacuum=3)
    atoms.get_potential_energy()
    # atoms.set_initial_magnetic_moments([0.5, 0.5])
    # calculator.set(charge=1)
    # atoms.get_potential_energy()

    # read again
    new_atoms = xtp.read_atoms('H2.orb')

    # assert
    assert(atoms.calc.results['energy'] == new_atoms.calc.results['energy'])