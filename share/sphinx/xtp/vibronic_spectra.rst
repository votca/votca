.. _vibronic-spectra:

Vibronic Spectra
################

When an electron is removed from a molecule (photoemission, UPS), added to
it (inverse photoemission, IPES) or excited (absorption), the nuclei respond
to the new charge distribution. A spectrum therefore does not consist of
sharp lines at the quasiparticle (QP) or excitation energies but of
vibrational progressions, whose shape and width are set by the coupling of
the electronic level to the molecular vibrations. VOTCA-XTP computes these
couplings from analytic GW and GW-BSE gradients and a ground-state Hessian,
and turns them into vibronically resolved spectra of single molecules and of
morphologies.

The tools are in the Python package ``pyxtp.vibronic`` and in the script
``xtp_vibronic``; both are installed with :doc:`pyxtp`. The tutorial
:doc:`../xtp-tutorials/VIBRONIC_SPECTRA` goes through all steps for
formaldehyde.


Overview
********

A vibronic spectrum needs three ingredients, all at the minimum of the
electronic ground state and with the same DFT settings:

* the **normal modes** (from the ground-state Hessian),
* the **gradients** of the energies whose spectrum is wanted (QP levels for
  UPS and IPES, BSE excitation energies for absorption), and
* the **vertical energies** at which the bands are placed: the QP or
  excitation energies of the molecule, or the site energies of a
  morphology (e.g. from QM/MM).

The work is split into steps that each leave a file behind. Only the first
three are expensive, and none of them has to be repeated to change
temperature, broadening, levels or plots:

============================== =============================== ============================================== ====================================================
Step                           Cost                            Result                                         Repeat when
============================== =============================== ============================================== ====================================================
1. ground-state geometry       one DFT optimization            ``minimum.xyz``                                the method changes
2. Hessian                     :math:`6N+1` DFT force runs     one directory per displacement, ``modes.npz``  never (interrupted runs resume)
3. GW(-BSE) with gradients     one GW(-BSE) run + gradients    orb file                                       never; further levels or states are added to it
4. spectrum                    seconds                         text file, one column per level                parameters, levels or site energies change
5. figure                      seconds                         PDF/SVG/PNG                                    broadening, axis range, experiment overlay change
============================== =============================== ============================================== ====================================================

On the command line:

.. code-block:: bash

   xtp_vibronic hessian -x minimum.xyz -w hessian -t 4 -o options.xml
   xtp_vibronic modes -w hessian -m modes.npz
   xtp_tools -e dftgwbse -o dftgwbse.xml -n 8   # with qp_gradients, see below
   xtp_vibronic spectrum -m modes.npz -g molecule.orb -p ionization -o ups.dat
   xtp_vibronic plot -s ups.dat --broadening 0.1 --onset -o ups.pdf

``xtp_vibronic <command> --help`` lists all options of a command.


Theory
******

Couplings
=========

The model is that of displaced harmonic oscillators: the initial state (the
neutral ground state) has the harmonic potential of its Hessian, and the
final state (cation, anion or excited state) has the same frequencies and
normal modes, displaced by the linear term of the energy gap
:math:`\Delta E = E_\mathrm{final} - E_\mathrm{initial}`. With the gradient
:math:`\mathbf g = \nabla \Delta E` at the ground-state minimum and the
mass-weighted normal modes :math:`\mathbf e_\nu` with frequencies
:math:`\omega_\nu`, the projections

.. math::

   q_\nu = \mathbf e_\nu \cdot \mathbf M^{-1/2} \mathbf g

give the reorganization energy and the Huang-Rhys factor of every mode,

.. math::

   \lambda_\nu = \frac{q_\nu^2}{2\omega_\nu^2}, \qquad
   S_\nu = \frac{\lambda_\nu}{\hbar\omega_\nu}.

:math:`\lambda = \sum_\nu \lambda_\nu` is the relaxation energy of the final
state from the vertical point: the vertical and the adiabatic (0-0) energy
differ by :math:`\lambda`. For charge transfer between two molecules the
reorganization energy is the sum of the :math:`\lambda` of charging and
discharging. The sign of the gradient does not matter.

For a QP level :math:`p`, :math:`\Delta E` is the QP energy
:math:`E_p` itself: removing an electron from level :math:`p` costs
:math:`-E_p`, adding one to an empty level gains it. For a BSE state it is
the excitation energy :math:`\Omega_n`. Both gradients are analytic and
include the response of the QP corrections (see
:doc:`excited_state_gradients`).

Translations and rotations are projected out of the Hessian exactly (Eckart
conditions at the reference geometry): there are always :math:`3N-6`
(:math:`3N-5` for linear molecules) vibrations, and no frequency threshold
is needed to remove the rigid-body motions.

Line shape
==========

The thermal line shape of the transition, with the initial state in thermal
equilibrium at temperature :math:`T`, is the Fourier transform of the
generating function (Lax, Kubo),

.. math::

   F(E) = \frac{1}{2\pi}\int d\tau\, e^{i(E-E_{00})\tau}
   \exp\Big(\sum_\nu S_\nu\big[(n_\nu+1)e^{-i\omega_\nu\tau}
   + n_\nu e^{i\omega_\nu\tau} - (2n_\nu+1)\big]\Big)
   \, e^{-i\lambda_\mathrm{cl}\tau - \lambda_\mathrm{cl}k_BT\tau^2}
   \, e^{-\sigma^2\tau^2/2 - \gamma|\tau|},

with :math:`\tau = t/\hbar`, :math:`\hbar = 1` in the exponents, and the
Bose occupations :math:`n_\nu`. At :math:`T=0` it is a Poisson progression
of lines :math:`E_{00} + k\hbar\omega_\nu` with weights
:math:`e^{-S}S^k/k!`; at high temperature it becomes the Marcus Gaussian.
The further terms are

* :math:`\lambda_\mathrm{cl}`: a classical reorganization energy, for
  modes that are not in the Hessian (the environment, *outer sphere*) or
  that are better treated classically, such as torsions that are hindered
  in a solid (option ``classical_below``). It adds a Gaussian of variance
  :math:`2\lambda_\mathrm{cl}k_BT`.
* :math:`\sigma` and :math:`\gamma`: Gaussian (standard deviation) and
  Lorentzian (half width) broadening of every vibronic line, e.g. for the
  lifetime of the final state.

:math:`F` is normalized, its mean is the vertical energy
:math:`E_\mathrm{vert} = E_{00} + \lambda`, and it is evaluated with one
fast Fourier transform. A warning is printed if the time integral has not
decayed on the grid (no broadening and no classical part).

Process and energy axis
=======================

Which side of :math:`E_{00}` the progression lies on depends on the
measurement and on the energy axis it is plotted on. ``pyxtp.vibronic``
names the physical process instead:

===================== ======================================================== ===========
process               energy axis (and site energies)                          progression
===================== ======================================================== ===========
``absorption``        photon energy, :math:`S_0\to S_n` (:math:`\Omega_n`)       above
``emission``          photon energy, :math:`S_n\to S_0` (:math:`\Omega_n`)       below
``electron-addition`` QP energy of the added electron (LUMO, LUMO+1, ...;      above
                      inverse photoemission)
``electron-removal``  QP energy of the removed electron (HOMO, HOMO-1, ...;    below
                      photoemission)
``ionization``        ionization energy, :math:`-E_\mathrm{QP}` of an           above
                      occupied level (photoemission, binding energy)
``electron-affinity`` electron affinity, :math:`-E_\mathrm{QP}` of an empty     below
                      level (inverse photoemission)
===================== ======================================================== ===========

For emission the modes should strictly be those of the excited state.

Ensembles and spectra of several levels
=======================================

A spectrum of a film is the sum of the line shapes of all molecules, each
placed at its own vertical energy :math:`E_i` (e.g. QP or excitation
energies from QM/MM):

.. math::

   I(E) = \sum_i w_i F(E - E_i),

with optional weights :math:`w_i`, e.g. :math:`e^{-z_i/\ell}` for the
surface sensitivity of photoemission (depth :math:`z_i`, inelastic mean free
path :math:`\ell`). The line shape :math:`F` is that of one molecule; the
disorder of the environment enters through the distribution of the
:math:`E_i`, the polarization of the environment through their mean.

For UPS and IPES the bands of all levels are added, each with its own
couplings. Often site energies are known only for the frontier levels;
the other levels then use the distribution of the nearest level that has
them, shifted rigidly by the difference of the molecular QP energies. All
levels count equally: photoemission cross sections are not included.

The onset of a spectrum is determined as in experiment: the tangent at the
steepest point of the outer edge of the first band is extrapolated to zero
intensity.


Step by step
************

1. Ground-state geometry
========================

The couplings are those at the ground-state minimum, so the geometry has to
be optimized tightly with the same functional and basis set as the
following steps, e.g. with analytic DFT forces through :doc:`pyxtp` and
ASE:

.. code-block:: python

   from ase.io import read, write
   from ase.optimize import BFGS
   from pyxtp.vibronic import xtp_ground_state_factory

   settings = {'dftpackage/functional': 'XC_HYB_GGA_XC_PBEH',   # PBE0
               'dftpackage/basisset': 'def2-tzvp',
               'dftpackage/auxbasisset': 'aux-def2-tzvp'}
   factory = xtp_ground_state_factory(settings, nthreads=8)
   atoms = read('molecule.xyz')
   atoms.calc = factory('.', None)
   BFGS(atoms).run(fmax=5e-4)
   write('minimum.xyz', atoms)

Functionals are given by their libxc names (``XC_HYB_GGA_XC_PBEH`` for
PBE0). The residual forces of the minimum are reported in step 2.

2. Hessian
==========

``FiniteDifferenceHessian`` differentiates analytic forces by central
differences: every atom is displaced along :math:`x`, :math:`y` and
:math:`z` in both directions (:math:`6N+1` calculations with the reference;
``scheme='forward'`` needs :math:`3N+1` but has a larger error). The
default step of 0.005 Å reproduces analytic Hessians to about
1 cm\ :sup:`-1` (see Validation).

.. code-block:: python

   from pyxtp.vibronic import FiniteDifferenceHessian

   fd = FiniteDifferenceHessian(read('minimum.xyz'), factory,
                                workdir='hessian', step=0.005)
   fd.run()                        # or fd.run(part=i, nparts=n)
   modes = fd.normal_modes()
   print(modes.summary())
   modes.save('modes.npz')

or

.. code-block:: bash

   xtp_vibronic hessian -x minimum.xyz -w hessian -t 8 -o options.xml
   xtp_vibronic status -w hessian
   xtp_vibronic modes -w hessian -m modes.npz

* Every displacement runs in its own directory of ``workdir`` and is
  claimed with a lock file. A repeated run skips what is done, so an
  interrupted run resumes, and several processes can work on the same
  ``workdir`` at the same time, on one machine or on cluster nodes that
  share the file system. ``--part i --nparts n`` (``run(part=i,
  nparts=n)``) splits the jobs deterministically, e.g. for array jobs. The
  reference geometry is computed first, and the displaced SCFs start from
  its orbitals.
* ``xtp_vibronic status --clear-locks`` removes the locks of jobs that
  crashed.
* The workdir remembers geometry, step, scheme and atoms; a different
  setup is refused instead of mixing results.
* ``xtp_vibronic modes`` reports the largest residual force at the
  reference geometry and the asymmetry of the finite-difference Hessian.
  Imaginary frequencies mean that the geometry is not a minimum; their
  modes are left out of the couplings with a warning.
* ``indices`` (``--atoms``) displaces only some atoms (partial Hessian).
* Any ASE calculator can be used through its own factory
  ``f(directory, displacement)``; an existing ORCA Hessian can be read with
  ``read_orca_hessian`` (``xtp_vibronic spectrum --orca-hess``).

3. Gradients
============

DFT and GW-BSE are run once at the minimum, with the gradients that are
needed, using ``dftgwbse`` options:

* ``qp_gradients``: QP levels for UPS and IPES, as a list of levels or
  ranges relative to HOMO and LUMO or as absolute MO indices, e.g.
  ``homo-9:lumo+9``, ``homo,lumo`` or ``30:40``. The levels must lie in
  the QP window (``gwbse.qpmin``, ``gwbse.qpmax``). No BSE is needed.
* ``excited_state_forces``: one BSE state for absorption or emission, e.g.
  ``s1`` (needs ``dftpackage.xtpdft.compute_forces``).

.. code-block:: python

   from pyxtp import xtp

   calc = xtp(nthreads=8, excited_state='s1')
   calc.set(**settings, **{'gwbse/gw/mode': 'evGW'})
   calc.options.qp_gradients = 'homo-9:lumo+9'
   atoms.calc = calc
   atoms.get_potential_energy()        # writes the orb file

The gradients are stored in the orb file. For more levels or another
excited state later, DFT and GW-BSE need not be repeated:

.. code-block:: bash

   xtp_vibronic gradients molecule.orb --levels homo-15:lumo+15 -o options.xml
   xtp_vibronic gradients molecule.orb --state s2 -o options.xml

This runs ``dftgwbse`` with ``tasks=gradients`` and must use the options of
the original run (basis sets and functional are compared with the file; a
different GW-BSE setup is caught by the QP and excitation energies rebuilt
for the gradient). Levels whose gradients are stored already are not
recomputed. From Python: ``pyxtp.vibronic.add_gradients``.

DFT, GW and the three-index integrals are shared by all levels; every
level then costs one coupled-perturbed Kohn-Sham (Z-vector) solution and
one pass over the integral derivatives. The
gradients are those of the perturbative QP energies; see
:doc:`excited_state_gradients` for the methods covered (closed shell, RI,
G0W0 or evGW with the plasmon-pole or exact self-energy, no embedding).
Gradients from other sources can be used as well: ``read_gap_gradient``
reads orb files, ORCA-style ``.engrad`` files and arrays in eV/Å.

4. Couplings and spectra
========================

Single molecule
---------------

.. code-block:: python

   from pyxtp.vibronic import NormalModes, QPLevels, level_spectrum

   modes = NormalModes.load('modes.npz')
   levels = QPLevels.from_orb('molecule.orb')
   for label, c in zip(levels.labels, levels.couplings(modes)):
       print(label, c.summary(3))

   ups = level_spectrum(modes, levels, 'ionization', temperature=300,
                        sigma=0.02)
   ups.save('ups.dat')

``level_spectrum`` uses all occupied levels with gradients for
``electron-removal`` and ``ionization`` and all empty ones for
``electron-addition`` and ``electron-affinity``; ``select`` restricts
them, e.g. ``select=['homo-2', 'homo-1', 'homo']``. Without site energies
every band is placed at the molecular QP energy. The command line
equivalent is

.. code-block:: bash

   xtp_vibronic spectrum -m modes.npz -g molecule.orb -p ionization \
       --sigma 0.02 -o ups.dat

For one transition, e.g. the absorption to S\ :sub:`1`:

.. code-block:: python

   from pyxtp.vibronic import (couplings_from_gradient, read_gap_gradient,
                               vibronic_lineshape)

   s1 = couplings_from_gradient(modes, read_gap_gradient('molecule.orb'))
   shape = vibronic_lineshape(s1, 300, 'absorption', sigma=0.02)
   # shape.energy is relative to the vertical excitation energy,
   # shape.zero_zero is E_00 - E_vert

(``xtp_vibronic lineshape -g molecule.orb -p absorption ...``).
``couplings_from_geometry`` takes an optimized excited-state geometry
instead of the gradient; the difference between the two routes measures
anharmonicity and frequency changes.

Morphology
----------

Site energies are given per level, on the energy axis of the process:

.. code-block:: python

   film = level_spectrum(modes, levels, 'electron-addition', temperature=300,
                         lambda_classical=0.01,
                         site_energies={'lumo': e_lumo, 'lumo+1': e_lumo1},
                         weights=np.exp(-depth / 2.0))
   for label, note in film.notes.items():
       print(label, note)

.. code-block:: bash

   xtp_vibronic spectrum -m modes.npz -g molecule.orb -p electron-addition \
       --lambda-classical 0.01 --sites sites.dat --columns lumo=3,lumo+1=4 \
       --where '1>4.0' --weight-column 7 -o ipes.dat

``--columns`` assigns columns (0-based) of the site file to levels,
``--where`` selects rows (e.g. sites more than 4 nm below the surface, if
column 1 holds the depth). Levels without site energies are shifted copies
of the nearest level that has them; ``notes`` (and the header of the
output file) say how every level was placed.

Parameters
----------

======================== ==================================================================
``temperature``          K; enters the Bose occupations and the classical term
``lambda_classical``     classical (outer-sphere) reorganization energy, eV
``sigma``, ``gamma``     Gaussian (standard deviation) and Lorentzian (half width)
                         broadening of every vibronic line, eV
``classical_below``      modes below this energy are moved to the classical term
                         (``--classical-below`` in cm\ :sup:`-1`)
``min_wavenumber``       modes below this frequency (cm\ :sup:`-1`) are left out;
                         their :math:`\lambda` is reported
``step``                 energy grid, eV
======================== ==================================================================

The instrument resolution is best applied when plotting (below), so that
it can be changed without recomputing the spectrum.

5. Figures
==========

``pyxtp.vibronic.plot`` (needs ``matplotlib``) makes figures for
publications: vector output, journal column widths (``SINGLE_COLUMN``
3.375 in, ``DOUBLE_COLUMN`` 7 in), and colours that are distinguishable
with colour-vision deficiencies. Every level keeps its colour in all
figures (HOMO and LUMO first, then by the distance from the gap).

.. code-block:: python

   from pyxtp.vibronic import Spectrum
   from pyxtp.vibronic.plot import paper_style, plot_spectrum, plot_couplings

   with paper_style():
       fig = plot_spectrum(Spectrum.load('ipes.dat'), broadening=0.3,
                           onset=True, experiment=(e_exp, i_exp))
       fig.savefig('ipes.pdf')
       print(fig.onset)                  # onset and first peak, eV
       plot_couplings(levels.couplings(modes)[0], title='LUMO').savefig('lumo.pdf')

.. code-block:: bash

   xtp_vibronic plot -s ipes.dat --broadening 0.3 --onset \
       --experiment measured.dat -o ipes.pdf
   xtp_vibronic lineshape -m modes.npz -g molecule.orb -p absorption \
       --couplings s1_modes.txt -o s1.dat
   xtp_vibronic plot -c s1_modes.txt --title S1 -o s1_modes.pdf

``plot_spectrum`` draws the contributions of the levels, the intrinsic
total, the total with the instrument broadening (Gaussian, standard
deviation) and an optional measured spectrum, scaled to the broadened
total. With ``onset`` the tangent construction and the first peak are
marked. ``plot_couplings`` shows :math:`\lambda_\nu` (or :math:`S_\nu`)
per mode and the running sum.


Validation
**********

* Hessian: water, PBE0/def2-SVP, finite differences of xtp forces
  1642/3866/3966 cm\ :sup:`-1`, analytic Hessian of PySCF
  1643.5/3865.2/3964.9 cm\ :sup:`-1`.
* QP gradients against finite differences of the QP energies with GW
  recomputed at displaced geometries (water, HOMO-2 to LUMO+2): relative
  errors below :math:`10^{-6}` for G0W0 and evGW with the plasmon-pole and
  the exact self-energy.
* Line shapes: Poisson progression at :math:`T=0`, Marcus limit at high
  temperature, normalization, mean and variance, detailed balance
  (unit tests in ``test_vibronic.py``).


Limitations and good practice
*****************************

* **Linear coupling only.** The final state has the frequencies and modes
  of the initial state: no frequency changes, no Duschinsky rotation, and
  modes that are not displaced by symmetry do not appear even if their
  frequency changes (e.g. the out-of-plane mode of formaldehyde in
  S\ :sub:`1`). Herzberg-Teller effects are not included.
* **Consistency.** Geometry, Hessian and gradients must come from the same
  functional and basis set, and the geometry must be a tight minimum: the
  gradient of the gap is only that of the final state if the ground-state
  gradient vanishes.
* **Soft modes.** Low-frequency torsions have large Huang-Rhys factors but
  small :math:`\lambda`; in a solid they are hindered. ``classical_below``
  (e.g. 100 cm\ :sup:`-1`) moves them into the classical term,
  ``min_wavenumber`` leaves them out.
* **Intensities.** All levels and sites count equally; cross sections of
  photoemission and the oscillator strengths of other transitions are not
  included.
* **QP levels.** The gradients are those of the perturbative QP energies
  (diagonal QP Hamiltonian). Levels with two QP solutions of similar
  weight (e.g. high virtual levels in evGW) should be kept out of the QP
  window. For degenerate levels the gradient of their mean energy is used.
