#!/usr/bin/env python3
"""Example to optimize the geometry of CO"""
from pyxtp import xtp
from ase.io import write
from ase.build import molecule
from ase.optimize import QuasiNewton

atoms = molecule('CO')
atoms.rattle()

calc = xtp(nthreads=2)
calc.select_force(energy='singlets', level=0, dynamic=False)
atoms.calc = calc

atoms.get_forces()

dyn = QuasiNewton(atoms, trajectory='test.traj')
dyn.run(fmax=0.01)
write('final.xyz', atoms)
