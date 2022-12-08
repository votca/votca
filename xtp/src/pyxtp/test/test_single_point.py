def test_single_point():
    import unittest

    from ase.build import molecule
    try:
        from pyxtp.ase_calculator.ase_calculator import xtp 
    except ImportError:
        raise unittest.SkipTest('xtp not available')

    
    if 1:
        calculator = xtp()
        atoms = molecule('H2', calculator=calculator)
        atoms.center(vacuum=3)
        atoms.get_potential_energy()
        # atoms.set_initial_magnetic_moments([0.5, 0.5])
        # calculator.set(charge=1)
        # atoms.get_potential_energy()

    # read again
    # t = io.read(txt, index=':')
    # assert isinstance(t, list)
    # M = t[1].get_magnetic_moments()
    # assert abs(M - 0.2).max() < 0.1
