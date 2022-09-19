from pyxtp import DFTGWBSE, Molecule, Visualization

mol = Molecule()
mol.add_atom("C", 0, 0, 0)
mol.add_atom("O", 1.2, 0, 0)

dft = DFTGWBSE(mol)
dft.options.data.dftpackage.basisset = 'dzp'
dft.run()