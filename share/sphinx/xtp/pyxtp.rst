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

Analytic excited-state forces
*****************************

Theory, the dftgwbse options, recommended settings and validation are in
:doc:`excited_state_gradients`.

``xtp(excited_state='s1')`` (or ``'t1'``, ...) makes ``get_potential_energy``
return the total energy of that BSE state (DFT total energy plus the static
excitation energy) and ``get_forces`` its analytic forces, computed after the
BSE by ``dftgwbse.excited_state_forces`` (``dftpackage.xtpdft.compute_forces``
is switched on automatically). This needs one GW-BSE run per geometry instead
of one per displaced coordinate:

.. code-block:: python

   from ase.optimize import BFGS
   calc = xtp(nthreads=4, excited_state='s1')
   calc.set(**{'dftpackage/basisset': 'def2-svp',
               'dftpackage/auxbasisset': 'aux-def2-svp'})
   atoms.calc = calc
   BFGS(atoms).run(fmax=0.01)

The gradient (``DFTEngine::BSEExcitationGradient``) is that of the BSE
excitation energy (full or TDA) with static RPA screening, including the
response of the Kohn-Sham orbitals, of the screening and of all integrals.
With the plasmon-pole model (``gw.sigma_integrator`` ppm, the default) and
the full-frequency self-energy (exact) it also includes the response of the
QP corrections, for G0W0 and for evGW
(an adjoint equation through the self-consistent energies, including
``gw.out_of_window_shift``), so the forces are the derivative of the energy
returned (checked against finite differences to about 1e-4). With the
exact self-energy, QP solutions can jump between branches near poles of the
self-energy, which makes the energy discontinuous along a path; the forces
belong to the branch found at that geometry. For evGW,
energy and forces agree to the convergence of the GW iterations
(``gw.sc_limit``). For geometry optimizations with evGW, set
``gw.out_of_window_shift`` to ``nearest``: the default ``max`` shifts the
levels outside the GW window by the largest correction in the window, so
the energy has kinks where that level changes, and optimizations can stall
there (formaldehyde S1, PBE0/def2-TZVPP: no convergence in 80 steps with
``max``, 8 steps with ``nearest``). For contour deformation (cda)
the QP corrections are kept fixed (cf. Holzer, arXiv:2611.12617); for water
(S1, PBE0/def2-SVP) this neglects about 10 % of the forces, and the energy,
which contains the full GW response, is then not exactly consistent with
the forces, so use an optimizer driven by the forces (e.g. ASE's ``BFGS``
or ``FIRE``). Requirements: xtp DFT
with an auxiliary basis, closed shell, ``bse.use_Hqp_offdiag`` false, no QSGW, no embedding or environment
screening. The run stops with an error if the excitation energy rebuilt for
the gradient does not reproduce the BSE eigenvalue.

Root following: with ``xtp(excited_state='s2', follow_root=True)`` the state
is followed from geometry to geometry instead of keeping its label. From the
second geometry on, ``dftgwbse.excited_state_forces_follow`` points to the orb
file of the previous geometry, and energy and forces belong to the state of
the same spin with the largest exciton overlap
:math:`|\langle X|X'\rangle - \langle Y|Y'\rangle|` with the previous state
(excitons in the AO basis, with the overlap of the basis functions at the two
geometries); ``calc.tracked_state`` gives its current label. The overlaps
are logged; a warning flags steps where the largest one is below 0.5 or the
next largest is above 80 % of it. Example, formaldehyde (G0W0/PBE0,
def2-SVP) along the C=O stretch: the bright state that is S4 at 1.20 Å is
followed through crossings to S3 (1.35 Å) and S2 (1.45 Å), with overlaps
above 0.97 for steps of 0.05 Å.

Excited-state dipoles by finite fields
**************************************

``pyxtp.finite_field.FiniteFieldDipoles`` computes the change of the dipole
moment on excitation, :math:`\Delta\mu_n = \mathrm{d}\Omega_n/\mathrm{d}F`,
and the diagonal of :math:`\Delta\alpha_n = -\mathrm{d}^2\Omega_n/\mathrm{d}F^2`
from GW-BSE excitation energies :math:`\Omega_n` in small uniform fields
(``dftpackage.externalfield``), together with the ground-state dipole and
polarizability from the DFT energy. Unlike the dipole of the excited-state
density (``Orbitals::DensityMatrixFull``, which is the unrelaxed difference
density :math:`XX^T + YY^T`), these are *relaxed*: they include the response
of the orbitals, of the quasiparticle energies and of the screened
interaction. For charge-transfer states the two can differ by a factor of
about two.

.. code-block:: python

   from ase.io import read
   from pyxtp import xtp
   from pyxtp.finite_field import FiniteFieldDipoles, write_xml

   calc = xtp()
   calc.set(basis='def2-tzvp', **{'dftpackage/auxbasisset': 'aux-def2-tzvp'})
   ff = FiniteFieldDipoles(read('molecule.xyz'), options=calc.options,
                           nthreads=4, components='xyz', nstates=5)
   result = ff.run()
   print(result.table())
   write_xml(result, "dipoles.xml")

* Each field component needs two ``dftgwbse`` runs (``order=2``, field step
  ``field=1e-3`` Hartree/bohr) or four (``order=4``, which also reports a
  linearity check), plus one run without field. The runs are kept in
  ``workdir`` and reused if the calculation is repeated, e.g. with more
  components.
* The SCF is converged with ``convergence_tightness=verytight``; the field
  changes the energies by only about :math:`10^{-4}` Hartree. Self-consistent
  GW (evGW) must be well converged, too: a GW window whose highest levels
  have two quasiparticle solutions can make the differences noisy (use an
  explicit ``qpmax``/``bsemax`` below such levels).
* States are matched between the runs by the overlap of their BSE amplitudes,
  after transforming them into the zero-field orbitals, so that states which
  cross in the field are followed correctly. States that mix with others, or
  that are matched with an overlap below ``overlap_threshold``, are marked;
  for groups of mixing states the derivative of the sum of their energies is
  given as well, which is reliable even when the individual values are not.
* Conventions: atomic units; ``externalfield`` adds :math:`-\mathbf{F}\cdot\mathbf{r}`
  to the one-electron Hamiltonian, so that :math:`\mathrm{d}E/\mathrm{d}\mathbf{F}`
  is the dipole moment. The ground-state dipole is given with respect to the
  centre of mass (the origin of the field operator).

The command-line script ``xtp_excited_dipoles`` wraps this class. It is
installed only together with ``pyxtp`` (see Installation) and has to run with
the Python interpreter ``pyxtp`` was installed for. Example:
``xtp_excited_dipoles -x molecule.xyz -o options.xml -t 4 --xml result.xml``,
where ``options.xml`` holds ``dftgwbse`` options as for ``xtp_tools``.

Vibronic spectra
****************

``pyxtp.vibronic`` and the script ``xtp_vibronic`` compute vibronically
resolved photoemission, inverse photoemission and absorption spectra from a
finite-difference Hessian of analytic DFT forces and analytic GW and GW-BSE
gradients, for single molecules and morphologies, and make figures for
publications. See :doc:`vibronic_spectra` and the tutorial
:doc:`../xtp-tutorials/VIBRONIC_SPECTRA`.

Examples
********

See the tutorials :doc:`../xtp-tutorials/single_point` (energy and absorption
spectrum of CO) and :doc:`../xtp-tutorials/geometry_optimization` (excited-state
forces and geometry optimization of CO) and
:doc:`../xtp-tutorials/VIBRONIC_SPECTRA` (vibronic UPS, IPES and absorption
spectra of formaldehyde).
