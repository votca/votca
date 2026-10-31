XTP ASE Interface
#####################

The Python package ``pyxtp`` provides an `ASE <https://wiki.fysik.dtu.dk/ase/>`__
calculator ``xtp`` that runs DFT and GW-BSE calculations (the ``dftgwbse``
tool) on an ``ase.Atoms`` object.

Installation
**************

The interface is built with VOTCA-XTP if ``pybind11`` and the Python
development files are found by CMake. ``make install`` installs the compiled
bindings and, if the Python modules ``ase``, ``xmltodict``, ``h5py``,
``lxml`` and ``rdkit`` are available, the ``pyxtp`` package into the Python
user site-packages directory (set ``CMAKE_INSTALL_PYTHON_LIBDIR`` to choose
another directory). The plotting functions of ``Visualization`` also need
``matplotlib``.

Environment variables
=====================

The environment variable :envvar:`VOTCASHARE` must be defined in order for the
interface to find the default xml files and create the options data structure.
You can set this environment variable in your shell configuration file
(``VOTCARC.bash`` sets it, too):

.. highlight:: bash

::

  $ export VOTCASHARE=<installation prefix>/share/votca

.. highlight:: python

Or within python itself:

  >>> os.environ['VOTCASHARE'] = '<installation prefix>/share/votca'

Options
*****************

All options of the ``dftgwbse`` tool are available in the ``options``
attribute of the calculator (an instance of ``pyxtp.options.Options``),
which mirrors the XML option files, e.g. ``calc.options.dftpackage.basisset``
or ``calc.options.gwbse.gw.mode``. Options that are not changed keep their
default values (see :doc:`dftgwbse`). XC functionals are given by their
libxc names, e.g. ``XC_HYB_GGA_XC_PBEH`` (the default) or
``XC_GGA_X_PBE XC_GGA_C_PBE``.

XTP calculator
*****************

``xtp(label=None, nthreads=1, options=None, directory='./', **kwargs)``

-  ``nthreads``: number of threads of the calculation
-  ``options``: an ``Options`` object to start from (default: the defaults
   of ``dftgwbse``)
-  further keyword arguments are passed to ``set``

Frequently used options can be specified with the calculator's ``set``
routine (or as keyword arguments of the constructor):

==================== ========= ============================= =================================
keyword              type      default value                 option
==================== ========= ============================= =================================
``xc``               ``str``   ``XC_HYB_GGA_XC_PBEH``        ``dftpackage.functional``
``basis``            ``str``   ``def2-tzvp``                 ``dftpackage.basisset``
``charge``           ``int``   0                             ``dftpackage.charge``
==================== ========= ============================= =================================

Other options can be set with ``set`` by their path, e.g.
``calc.set(**{'dftpackage/auxbasisset': 'aux-def2-svp'})``.

Results
*******

Besides ``get_potential_energy()`` and ``get_forces()`` (ASE), the calculator
provides e.g. ``get_dft_energy()``, ``get_qp_corrections()`` and
``get_oscillator_strengths(dynamic=False)``. ``pyxtp.Visualization(atoms)``
plots quasiparticle corrections, densities of states and absorption spectra.

Forces
******

``select_force(energy='energy', level=0, dynamic=False)`` chooses the energy
whose forces ``get_forces`` (and therefore an ASE optimizer) uses:

==================== ===============================================================
``energy``           energy
==================== ===============================================================
``'energy'``         DFT total energy (analytic forces if the DFT engine provides
                     them, otherwise finite differences)
``'ks'``             total energy with Kohn-Sham level ``level`` excited
``'qp'``             total energy with quasiparticle level ``level`` excited
``'singlets'``       BSE singlet ``level`` (``dynamic=True``: with dynamical
                     screening)
``'triplets'``       BSE triplet ``level`` (``dynamic`` as for singlets)
==================== ===============================================================

All energies except ``'energy'`` use finite differences of the total energy
(displacement ``eps=0.001`` Å in ``get_forces``), i.e. one full DFT and
GW-BSE calculation per displaced coordinate.

Examples
********

See the tutorials :doc:`../xtp-tutorials/single_point` (energy and absorption
spectrum of CO) and :doc:`../xtp-tutorials/geometry_optimization` (excited-state
forces and geometry optimization of CO).
