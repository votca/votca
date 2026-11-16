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

Vibronic line shapes
********************

``pyxtp.vibronic`` turns a ground-state Hessian and the gradient of an energy
gap into Huang-Rhys factors and thermal vibronic line shapes, and sums them
over the vertical energies of a morphology (e.g. QP levels or excitation
energies from QM/MM) to give ensemble spectra.

Model: displaced harmonic oscillators with the frequencies of the initial
state. With :math:`g` the gradient of :math:`E_\mathrm{final}-E_\mathrm{initial}`
at the initial-state minimum, the mass-weighted normal-mode components
:math:`q_\nu = e_\nu \cdot M^{-1/2} g` give

.. math::

   \lambda_\nu = \frac{q_\nu^2}{2\omega_\nu^2}, \qquad
   S_\nu = \frac{\lambda_\nu}{\hbar\omega_\nu},

and the line shape is the Fourier transform of the generating function

.. math::

   F(E) = \frac{1}{2\pi}\int d\tau\, e^{i(E-E_{00})\tau}
   \exp\Big(\sum_\nu S_\nu\big[(n_\nu+1)e^{-i\omega_\nu\tau}
   + n_\nu e^{i\omega_\nu\tau} - (2n_\nu+1)\big]\Big)
   e^{-i\lambda_\mathrm{cl}\tau - \lambda_\mathrm{cl}k_BT\tau^2}
   e^{-\sigma^2\tau^2/2 - \gamma|\tau|},

with :math:`\tau=t/\hbar`, Bose occupations :math:`n_\nu`, a classical
(outer-sphere or soft-mode) reorganization energy :math:`\lambda_\mathrm{cl}`
and optional Gaussian (:math:`\sigma`) and Lorentzian (:math:`\gamma`)
broadening. It is evaluated with one FFT, is normalized, and is stored
relative to the vertical energy :math:`E_\mathrm{vert}=E_{00}+\lambda`,
which is its mean.

.. code-block:: python

   from ase.io import read
   from pyxtp.vibronic import (FiniteDifferenceHessian, xtp_ground_state_factory,
                               couplings_from_gradient, read_gap_gradient,
                               vibronic_lineshape, ensemble_spectrum)

   settings = {'dftpackage/functional': 'XC_HYB_GGA_XC_PBEH',  # PBE0
               'dftpackage/basisset': 'def2-svp',
               'dftpackage/auxbasisset': 'aux-def2-svp'}
   fd = FiniteDifferenceHessian(read('minimum.xyz'),
                                xtp_ground_state_factory(settings, nthreads=4),
                                workdir='hessian', step=0.005)
   fd.run()                       # 6N + 1 force calculations, resumable
   modes = fd.normal_modes()      # translations/rotations projected out
   modes.save('modes.npz')

   # gradient of Omega_S1 from a dftgwbse run with excited_state_forces=s1
   c = couplings_from_gradient(modes, read_gap_gradient('s1.orb'),
                               min_wavenumber=10)
   print(c.summary())
   shape = vibronic_lineshape(c, temperature=300, process='absorption',
                              lambda_classical=0.05)
   energies, spectrum = ensemble_spectrum(site_energies, shape)

* **Hessian.** ``FiniteDifferenceHessian`` differentiates analytic forces
  (central differences, 6N + 1 calculations; ``scheme='forward'`` 3N + 1).
  Every displacement runs in its own directory of ``workdir`` and is claimed
  with a lock file: interrupted runs resume, several processes can work on
  the same ``workdir``, and ``run(part=i, nparts=n)`` splits the jobs for
  array jobs. ``xtp_ground_state_factory`` runs only the DFT part of
  ``dftgwbse`` with analytic forces and ``verytight`` SCF, and starts the
  displaced SCFs from the converged orbitals of the reference geometry. Any
  other ASE calculator works as well. ``indices`` restricts the
  displacements to some atoms (partial Hessian). The reference geometry
  should be a tightly converged minimum of the same method; the residual
  forces and the asymmetry of the finite-difference Hessian are reported by
  ``xtp_vibronic modes``. Existing Hessians can be read with
  ``read_orca_hessian``.
* **Gradients.** ``read_gap_gradient`` reads the gradient of
  :math:`\Omega_n` from an orb file with ``excited_state_forces`` (see
  above; :math:`F_\mathrm{ground}-F_\mathrm{excited}`), ORCA-style
  ``.engrad`` files or arrays in eV/Å. ``couplings_from_forces`` takes the
  forces of both states; ``couplings_from_geometry`` uses an optimized
  final-state geometry instead (the difference to the gradient route
  measures anharmonicity and frequency changes).
* **Process.** ``process`` (``--process``) sets the energy axis of the line
  shape and of the site energies, and with it the side of :math:`E_{00}` on
  which the vibronic progression lies:

  ===================== ================================================ ===========
  process               energy axis (site energies)                      progression
  ===================== ================================================ ===========
  ``absorption``        photon energy :math:`S_0\to S_n` (:math:`\Omega_n`) above
  ``emission``          photon energy :math:`S_n\to S_0` (:math:`\Omega_n`) below
  ``electron-addition`` QP level of the added electron (LUMO, ...;       above
                        inverse photoemission)
  ``electron-removal``  QP level of the removed electron (HOMO, ...;     below
                        photoemission)
  ``ionization``        ionization energy :math:`-E_\mathrm{HOMO}`        above
  ``electron-affinity`` electron affinity :math:`-E_\mathrm{LUMO}`        below
  ===================== ================================================ ===========

  The gradient is that of :math:`\Omega_n` or of the QP level involved; its
  sign does not matter (:math:`\lambda_\nu \propto q_\nu^2`). For emission
  the modes should strictly be those of the excited state.
* **Soft modes.** Modes below ``min_wavenumber`` are left out (their
  :math:`\lambda` is reported); ``classical_below`` treats soft modes as part
  of the classical term instead, which is appropriate for torsions hindered
  in a solid. Imaginary modes are always left out with a warning.
* **Limitations.** Linear coupling only: the final state has the frequencies
  and modes of the initial state (no frequency changes, no Duschinsky
  rotation), so modes that are not displaced by symmetry do not appear even
  if their frequency changes (e.g. the out-of-plane mode of formaldehyde in
  S1). Herzberg-Teller (intensity borrowing) effects are not included.

Spectra of several levels (UPS, IPES)
=====================================

With ``qp_gradients`` (e.g. ``homo-5:lumo+5``) ``dftgwbse`` stores the
analytic gradients of the QP energies of a range of levels (see
:doc:`excited_state_gradients`). ``QPLevels.from_orb`` reads them, and
``level_spectrum`` builds the whole spectrum: for every level that fits
the process -- occupied levels for ``electron-removal`` and ``ionization``
(photoemission), empty levels for ``electron-addition`` and
``electron-affinity`` (inverse photoemission) -- it computes the couplings
and the line shape and places it at the vertical energies of that level.

.. code-block:: python

   from pyxtp.vibronic import NormalModes, QPLevels, level_spectrum
   from pyxtp.vibronic.plot import paper_style, plot_spectrum

   modes = NormalModes.load('modes.npz')
   levels = QPLevels.from_orb('molecule.orb')     # dftgwbse with qp_gradients
   spectrum = level_spectrum(modes, levels, 'electron-addition',
                             temperature=300, lambda_classical=0.01,
                             site_energies={'lumo': e_lumo, 'lumo+1': e_lumo1})
   spectrum.save('ipes.dat')
   with paper_style():
       fig = plot_spectrum(spectrum, broadening=0.3, onset=True,
                           experiment=(e_exp, i_exp))
       fig.savefig('ipes.pdf')

* Without ``site_energies`` the spectrum is that of one molecule at its QP
  energies. With site energies of a morphology (per level, on the energy
  axis of the process), each level is summed over its sites. Levels
  without site energies use the distribution of the nearest level that has
  them, shifted rigidly by the difference of the molecular QP energies;
  ``Spectrum.notes`` records how every level was placed.
* All levels count equally; photoemission cross sections are not included.
* More levels (or the forces of another excited state) can be added later
  to the same orb file without repeating DFT and GW-BSE:
  ``add_gradients('molecule.orb', qp_levels='homo-9:lumo+9',
  options=calc.options)`` runs ``dftgwbse`` with ``tasks=gradients``;
  levels that are already stored are not recomputed.
* ``Spectrum.broadened(sigma)`` applies a Gaussian instrument broadening,
  ``Spectrum.onset()`` the tangent construction for the onset.
* ``pyxtp.vibronic.plot`` (needs ``matplotlib``) makes figures for
  publication: ``plot_spectrum`` (contributions of the levels, total,
  broadened total, measured spectrum, onset and first peak) and
  ``plot_couplings`` (:math:`\lambda_\nu` or :math:`S_\nu` per mode and
  their running sum), in ``paper_style()`` (journal column widths, vector
  fonts, colour-blind-safe colours in a fixed order).

The script ``xtp_vibronic`` provides the same steps on the command line
(installed together with ``pyxtp``):

.. code-block:: bash

   xtp_vibronic hessian -x minimum.xyz -w hessian -t 4 \
       --set dftpackage.functional=XC_HYB_GGA_XC_PBEH --set dftpackage.basisset=def2-svp \
       --set dftpackage.auxbasisset=aux-def2-svp     # repeatable, or -o options.xml
   xtp_vibronic hessian -w hessian --part 3 --nparts 16  # e.g. in an array job
   xtp_vibronic status -w hessian
   xtp_vibronic modes -w hessian -m modes.npz
   xtp_vibronic lineshape -m modes.npz -g s1.orb -p absorption -T 300 --lambda-classical 0.05 \
       --sites site_energies.dat --column 2 --where '1>4.0' -o spectrum.dat
   # all empty levels with QP gradients; LUMO and LUMO+1 from the morphology,
   # higher levels as shifted copies
   xtp_vibronic spectrum -m modes.npz -g molecule.orb -p electron-addition \
       --lambda-classical 0.01 --sites site_energies.dat \
       --columns lumo=3,lumo+1=4 --where '1>4.0' -o ipes.dat
   # later: more levels from the same orb file, no new DFT or GW
   xtp_vibronic gradients molecule.orb --levels homo-9:lumo+9 -o options.xml -t 8
   xtp_vibronic plot -s ipes.dat --broadening 0.3 --onset \
       --experiment measured.dat -o ipes.pdf
   xtp_vibronic plot -c s1_modes.txt --title S1 -o s1_modes.pdf

Examples
********

See the tutorials :doc:`../xtp-tutorials/single_point` (energy and absorption
spectrum of CO) and :doc:`../xtp-tutorials/geometry_optimization` (excited-state
forces and geometry optimization of CO).
