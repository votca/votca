#!/usr/bin/env python
"""Example to compute energies using XTP."""
from pyxtp import xtp, Visualization
from ase import Atoms


from pyxtp.xml_editor import (
    remove_empty_text_elements,
    insert_linked_xmlfiles,
    update_xml_text,
)

"""Run energy workflow."""
# define a molecule
atoms = Atoms('CO', positions=([0,0,0],[1.4,0,0]))

# define the calculator
calc = xtp(nthreads=2)

# change basis sets to a smaller one
# calc.options.dftpackage.basisset = 'def2-svp'
calc.options.dftpackage.basisset = 'XXX'
calc.options.dftpackage.auxbasisset = 'aux-def2-svp'

options = calc.options._opts.to_flat_dict()

# update the options from the user provided input
new_tree = update_xml_text(calc.options._xml_tree, options)

# remove the empty elements
new_tree_2 = remove_empty_text_elements(new_tree)

# write the file
new_tree_2.write('test.xml')

# calc.options._write_xml('test.xml')