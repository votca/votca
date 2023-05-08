XTP ASE Interface
#####################

Installation
**************


Environment variables
=====================

The environment variable :envvar:`VOTCASHARE` must be defined in order for the 
interface to find the default xml files and create the options data structure.
You can set this environment variable in your shell configuration file:

.. highlight:: bash

::

  $ export VOTCASHARE=/share/votca/

.. highlight:: python

Or within python itself:

  >>> os.environ['VOTCASHARE'] = '/share/votca/'

Option parser
*****************

The `XTPOptions` class is the preferred handler for setting all the options of the xtp calculator.

XTP calculator 
*****************

Frequently used options can be specified using
the calculator's `set` routine.

==================== ========= ============= =====================================
keyword              type      default value description
==================== ========= ============= =====================================
``label``            ``str``   None          Name of input and output files
``xc``               ``str``   None          XC functional to be used
``basisset``         ``str``   None          Basis set to be used
``charge``             int     0             Charge of the molecule
==================== ========= ============= =====================================

Single Point Caclulation
***************************

Here is an example of setting up a calculation on a water molecule: ::

    # Set up a water molecule
    from ase.build import molecule
    wat = molecule('H2O')

    # Set up a XTP calculator
    from pyxtp import xtp
   
    
    # change options for the calculations
    calc = xtp(label='water')
    calc.options.dftpackage.functional = 'PBE'
    calc.options.dftpackage.basisset = 'def2-svp'
    calc.options.dftpackage.auxbasisset = 'aux-def2-svp'

    # attach the calculator
    wat.calc = calc

    # compute the energy terms
    wat.get_potential_energy()

.. highlight:: python


Geometry Optimization
***************************

Here is an example of setting up a geometry of a carbon monoxide 
molecule using forces calculated on the first singlet state: ::

    # Set up a water molecule
    from ase.build import molecule
    co = molecule('CO')
    co.rattle()

    # Set up a XTP calculator
    from pyxtp import xtp
   
    # change options for the calculations
    calc = xtp()
    calc.options.dftpackage.functional = 'PBE'
    calc.options.dftpackage.basisset = 'def2-svp'
    calc.options.dftpackage.auxbasisset = 'aux-def2-svp'

    # specify the wich forces we want to use
    calc.select_force(energy='singlets', level=0, dynamic=False)

    # attach the calculator
    co.calc = calc

    # compute the energy terms
    from ase.optimize import QuasiNewton
    from ase.io import write
    dyn = QuasiNewton(co, trajectory='test.traj')
    dyn.run(fmax=0.01)
    write('final.xyz', atoms)

.. highlight:: python
