.. _excited-state-gradients:

Excited-State Gradients and Geometry Optimization
#################################################

VOTCA-XTP computes analytic nuclear gradients of BSE excitation energies,
and from them the forces on the atoms in an excited state. One force
evaluation costs about two GW-BSE calculations, independent of the
number of atoms, instead of :math:`6N` for the finite differences used
in earlier excited-state optimizations with VOTCA-XTP [Caylak:2021]_. Together
with an optimizer (through :doc:`pyxtp` and ASE) this gives excited-state
geometries, Stokes shifts and adiabatic excitation energies with GW-BSE.

The gradient includes the response of the Kohn-Sham orbitals, of the
screening, of all integrals and, for G0W0 and evGW, of the quasiparticle
(QP) corrections themselves. The latter follows Holzer [Holzer:2026]_;
neglecting it shortens excited-state bonds systematically (by about
0.003–0.006 Å for the C–C bond of acetylene in its
:math:`^1\Sigma_u^-` state, and by up to 0.03 Å for the C=O bond of
formaldehyde in S\ :sub:`1`).


Theory
******

The total energy of BSE state :math:`n` is the ground-state DFT energy
plus the excitation energy,

.. math::

    E_n(\mathbf{R}) = E_\text{DFT}(\mathbf{R}) + \Omega_n(\mathbf{R}),

and its gradient is that of the ground state (available with
``dftpackage.xtpdft.compute_forces``) plus that of :math:`\Omega_n`.

Excitation energy and its partial derivatives
=============================================

In the resolution-of-identity (RI) representation the BSE Hamiltonian is
built from the three-index quantities
:math:`B^P_{pq} = \sum_R L_{RP} (R|pq)`, with :math:`(R|pq)` the
three-centre Coulomb integrals in the molecular orbital (MO) basis and
:math:`L = V^{-1/2}` the inverse square root of the auxiliary Coulomb
metric, and from the QP energies :math:`E_p` (BSE) and the energies
:math:`\varepsilon'_p` in the RPA screening. For the full BSE, with
amplitudes :math:`X, Y` normalized as :math:`X^TX - Y^TY = 1`,

.. math::

    \Omega = X^T A X + Y^T A Y + 2 X^T B Y,

so that by the Hellmann-Feynman theorem only the explicit dependence of
:math:`A` and :math:`B` on the geometry contributes:

.. math::

    \mathrm{d}\Omega = \sum_P g^P \cdot \mathrm{d}B^P
                       + \sum_p w_p\, \mathrm{d}E_p
                       + \sum_p w'_p\, \mathrm{d}\varepsilon'_p .

:math:`g^P = \partial\Omega/\partial B^P` collects the exchange, the
direct (screened) term and the response of the static RPA screening
:math:`\epsilon^{-1}` through the RPA transitions; :math:`w` and
:math:`w'` are the diagonal elements of the BSE difference density and of
the screening response density.

The change of :math:`B^P` has an integral part, the derivatives of the
three-centre integrals and of the metric at fixed MO coefficients (both
contracted directly, without storing derivative tensors), and an orbital
part from the response of the MO coefficients.

Orbital response: the Z-vector
==============================

The orbital part is handled as for TDDFT gradients
[Furche:2002]_. The quantity
:math:`Q_{rp} = \sum_{P,q} B^P_{rq} (g^P_{pq} + g^P_{qp})` determines the
rotations within the occupied and within the virtual space directly; the
occupied-virtual rotations follow from one coupled-perturbed Kohn-Sham
(Z-vector) equation

.. math::

    (A + B)_\text{KS}\, Z = L,

with the orbital Hessian of the ground state, i.e. the response of the
Coulomb, exact-exchange and exchange-correlation potentials. The relaxed
difference density :math:`\tilde{P}` and the energy-weighted density
:math:`\tilde{W}` then give

.. math::

    \frac{\mathrm{d}\Omega}{\mathrm{d}\mathbf{R}} =
    (\text{integral part})
    + \sum_{\mu\nu} \tilde{P}_{\mu\nu}
      \frac{\partial F_{\mu\nu}}{\partial \mathbf{R}}\Big|_D
    - \sum_{\mu\nu} \tilde{W}_{\mu\nu}
      \frac{\partial S_{\mu\nu}}{\partial \mathbf{R}} ,

where the Fock-matrix derivative at fixed density is split into the
one-electron part, RI-J, RI-K and the exchange-correlation part. The
latter is evaluated analytically on the grid, including the moving basis
functions and integration grid (second functional derivatives of the
functional; no third derivatives are needed because the BSE kernel
contains no exchange-correlation kernel).

Response of the QP corrections
==============================

The QP energies solve

.. math::

    E_p = \varepsilon_p + (1 - a_x)\Sigma^x_{pp} - v^{xc}_{pp}
          + \Sigma^c_{pp}(E_p; \mathbf{R}),

with :math:`a_x` the exact-exchange fraction of the functional. They
depend on the geometry through the KS energies, the integrals in
:math:`\Sigma^x`, :math:`v^{xc}` and :math:`\Sigma^c`, and, in evGW, through
the QP energies in the screening and in the poles of :math:`G`. The
weights :math:`w_p` of :math:`\mathrm{d}\Omega` on the QP energies are
turned into weights on these quantities by an adjoint equation:

* **G0W0**: :math:`v_p = w_p / (1 - \Sigma^{c\prime}_{pp}(E_p))`, the
  renormalization factor.
* **evGW**: a linear system of the size of the GW window,
  :math:`(1 - J)^T v = \bar{w}`, with :math:`J` the derivative of the
  self-energies with respect to the QP energies, through the poles of
  :math:`G` and through the screening. The levels outside the GW window
  enter through the out-of-window shift (``gw.out_of_window_shift``). The
  system is solved with GMRES, preconditioned with the part of :math:`J`
  that comes from the poles of :math:`G` alone.

A reverse pass of :math:`\sum_p v_p \Sigma^c_{pp}` then gives the
contributions to :math:`g^P`, to the KS energies and to the RPA energies:

* **Plasmon-pole model** (``sigma_integrator`` ``ppm``): through the pole
  sums, the PPM parameters, and the eigen-decomposition of
  :math:`\epsilon(0)` (eigenvalue and eigenvector derivatives) to
  :math:`\epsilon(0)` and :math:`\epsilon(i\nu)`.
* **Exact self-energy** (``sigma_integrator`` ``exact``): :math:`\Sigma^c`
  is written as a matrix function of the RPA matrix
  :math:`C = \Delta^{1/2}(\Delta + 4bb^T)\Delta^{1/2}` and differentiated
  with divided differences in its eigenbasis (Daleckii-Krein).

The contour-deformation integrator (``cda``) is not covered; with it the
QP corrections are kept fixed.

Gradients of quasiparticle energies
===================================

The same machinery gives the nuclear gradient of the energy of a single QP
level :math:`E_p`: the "excitation energy" is then :math:`E_p` itself, so
the partial derivatives :math:`g^P` of the BSE part vanish and the weight
on the QP energies is one for level :math:`p`. The adjoint of the QP
equation (with the evGW self-consistency, if used), the Z-vector and the
integral derivatives are unchanged. ``dftgwbse`` computes these gradients
with the option ``qp_gradients``, a list of levels or ranges such as
``homo-5:lumo+5``, ``homo,lumo`` or ``30:40`` (absolute, 0-based MO
indices); no BSE is needed. The three-index integrals and the self-energy
model are set up once, each level then needs its own Z-vector. The result
is stored in the orb file (``qp_gradients``, one block of
:math:`N_\mathrm{atoms}` rows per level in Hartree/bohr, and
``qp_gradient_levels``).

* The gradients are those of the perturbative QP energies (diagonal of the
  QP Hamiltonian); the levels must lie in the QP window.
* For KS-degenerate levels only the gradient of the mean energy is
  defined; each partner gets this gradient (logged).
* The restrictions of the excited-state forces apply (closed shell, RI,
  G0W0 or evGW with ``ppm`` or ``exact``; otherwise the gradient of the KS
  energy with fixed corrections is returned).

Gradients can be added to an existing orb file without repeating DFT and
GW-BSE, which matters for large molecules: ``dftgwbse`` with
``tasks=gradients`` (and the options of the run that wrote
``job_name.orb``) reads the file, computes what ``qp_gradients`` and/or
``excited_state_forces`` ask for and writes it back. QP levels whose
gradients are already stored are not recomputed, so a range can be
extended later; the forces of another excited state replace the stored
ones. Basis sets and functional are compared with the file, and the QP and
excitation energies rebuilt for the gradient reveal a different GW-BSE
setup. From Python: ``pyxtp.vibronic.add_gradients``; on the command line
``xtp_vibronic gradients``.

Validation against finite differences of the QP energies (GW recomputed
at displaced geometries; water, PBE0/def2-SVP, QP window 0-9, levels
HOMO-2 to LUMO+2): relative errors below :math:`10^{-6}` for G0W0 and
evGW with both self-energy models.

The main use is the vibronic structure of photoemission and inverse
photoemission spectra, see :doc:`vibronic_spectra`.

Root following
==============

Excited states change order along an optimization. With root following,
the state at a new geometry is the one of the same spin with the largest
overlap with the state of the previous geometry,

.. math::

    O_i = \left| \langle X_\text{ref} | X_i \rangle
               - \langle Y_\text{ref} | Y_i \rangle \right|,

with the excitons in their atomic-orbital representation and the overlap
of the basis functions at the two geometries (so the moving basis is
accounted for). :math:`O_i` drops for any displacement and is meant for
consecutive optimization steps.


Usage
*****

Relaxed densities
=================

The gradient contracts the derivatives of the one-electron Hamiltonian
with :math:`\tilde{P}` (the BSE difference density, the Z-vector and,
with the QP response, the response of the QP corrections). This is the
relaxed difference density: for any perturbation :math:`\lambda h'` of the
one-electron Hamiltonian at fixed basis (an electric field, point charges,
multipoles),

.. math::

    \frac{\mathrm{d}\Omega}{\mathrm{d}\lambda}
    = \sum_{\mu\nu} \tilde{P}_{\mu\nu} h'_{\mu\nu},

and the same holds for single QP levels. Unlike the unrelaxed density
:math:`D_0 + XX^T + YY^T` it contains the orbital relaxation, so its
dipole moment is the derivative of the state energy with respect to a
field. Computing it needs no derivative integrals.

* ``dftgwbse``: ``relaxed_density`` (e.g. ``s1``, ``t1``, or a QP level such
  as ``pqp5``) stores :math:`D_0 + \tilde{P}` in the orb file (dataset
  ``relaxed_density``, attribute ``relaxed_density_state``; for a QP level
  the density with the electron added or removed) and prints its dipole
  change from the ground state. The response of the QP corrections
  follows ``excited_state_forces_qp_response``. With ``tasks gradients`` it
  is computed from an earlier run.
* ``qmmm``: ``qmregion.relaxed_density true`` polarizes the polar regions
  with the relaxed density of the tracked state in the iterative scheme
  (not with ``environment_screening``).

For water (PBE0/def2-SVP, full BSE), the dipole change of the relaxed S1
density agrees with finite-field derivatives of the S1 energy to
:math:`4\times 10^{-7}` (G0W0) and :math:`3\times 10^{-6}` (evGW)
:math:`e a_0`; with the QP corrections kept fixed it differs by about 7%.

From Python (ASE)
=================

The usual way is the ASE calculator, see :doc:`pyxtp`:

.. code-block:: python

   from ase.io import read
   from pyxtp import xtp
   from sella import Sella

   atoms = read('molecule.xyz')
   calc = xtp(nthreads=8, excited_state='s1', follow_root=True)
   calc.set(**{'dftpackage/basisset': 'def2-tzvpp',
               'dftpackage/auxbasisset': 'aux-def2-tzvpp',
               'gwbse/gw/mode': 'evGW',
               'gwbse/gw/sigma_integrator': 'exact',
               'gwbse/gw/out_of_window_shift': 'nearest'})
   atoms.calc = calc
   Sella(atoms, order=0, internal=True).run(fmax=0.01)

``get_potential_energy`` returns the total energy of the state and
``get_forces`` its analytic forces; with ``follow_root=True``,
``calc.tracked_state`` is the label of the followed state at the current
geometry.

With dftgwbse
=============

The forces are computed at the end of a ``dftgwbse`` run with

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Option
     - Meaning
   * - ``excited_state_forces``
     - BSE state, e.g. ``s1`` or ``t1``. Needs
       ``dftpackage.xtpdft.compute_forces true`` (set automatically by
       pyxtp).
   * - ``excited_state_forces_qp_response``
     - Include the response of the QP corrections (default ``true``;
       G0W0 and evGW with ``ppm`` or ``exact``). ``false`` keeps them
       fixed.
   * - ``excited_state_forces_follow``
     - Orb file of an earlier geometry: compute the forces for the state
       with the largest overlap with the state stored there (root
       following).

The forces (Hartree/Bohr) and the state are stored in the orb file
(dataset ``excited_state_forces``, attribute
``excited_state_forces_state``) and printed in the log. The log also
reports three checks:

* the excitation energy rebuilt for the gradient against the BSE
  eigenvalue (a warning above 1e-6 Ha means a loosely converged Davidson
  solution; pyxtp sets ``bse.davidson.tolerance`` to ``strict``; above
  1e-4 Ha the run stops, since the GW-BSE setup is then not covered by the
  gradient);
* the QP energies rebuilt from the self-energy against the GW ones (stop
  above 1e-4 Ha, e.g. for an evGW run that has not converged);
* for evGW, the RPA energies rebuilt with the out-of-window shift.


Recommended settings
********************

* **evGW with the exact self-energy** for geometries. evGW solves the QP
  equations iteratively and gives smooth energy surfaces; G0W0 surfaces
  can be rough where a QP solution lies close to a pole of the
  self-energy, and the plasmon-pole model changes bond lengths by up to
  about 0.006 Å relative to the exact self-energy.
* **out_of_window_shift nearest** with evGW. The default ``max`` shifts
  the levels outside the GW window by the largest correction in the
  window; the energy then has kinks where that level changes, and
  optimizations can stall there (dftgwbse warns).
* **A GW window without problematic levels.** If evGW does not converge
  because a core or high virtual level has two QP solutions of similar
  weight (the log lists these), restrict the window with
  ``ranges explicit`` and ``qpmin``/``qpmax``; the gradient refuses QP
  energies that do not solve the QP equations.
* **An optimizer in internal coordinates**, e.g. Sella, rather than BFGS
  in Cartesian coordinates: Cartesian BFGS cannot break the symmetry of
  the starting geometry, and its fixed step size overshoots near state
  crossings.
* **Root following** whenever states are close at the starting geometry.
* **Stop near S**\ :sub:`1`/**S**\ :sub:`0` **crossings.** GW-BSE
  describes excitations as linear response from the ground state; when
  the excitation energy drops to about 1 eV or below, or the ground-state
  energy jumps between steps, the description breaks down and the
  geometry is no longer meaningful.


Validation
**********

The analytic gradients agree with finite differences of the full GW-BSE
energy (GW recomputed at every displaced geometry) to about
:math:`2\times10^{-5}` relative (water, PBE0/def2-SVP: G0W0 and evGW with
the plasmon-pole model, evGW with the exact self-energy).

S\ :sub:`1` geometries of ten molecules of the 0–0 benchmark set of Loos
and Jacquemin [Loos:2019]_, optimized from the (EOM-)CCSD/def2-TZVPP
excited-state minima (PBE0, def2-TZVPP, full BSE, exact self-energy, GW
window of all occupied and the lowest ten virtual orbitals). Mean
absolute deviation of the heavy-atom bond lengths from the reference, in
Å, compared with the Turbomole implementation of Holzer [Holzer:2026]_:

.. list-table::
   :header-rows: 1

   * - State
     - evGW (XTP)
     - evGW (Holzer)
     - G0W0 (XTP)
     - G0W0 (Holzer)
   * - formaldehyde
     - 0.013
     - 0.014
     - 0.006
     - 0.001
   * - acetaldehyde
     - 0.007
     - 0.006
     - 0.013
     - 0.018
   * - acetylene
     - 0.009
     - 0.008
     - 0.005
     - 0.002
   * - HCN
     - 0.006
     - 0.006
     - 0.002
     - 0.006
   * - thioformaldehyde
     - 0.007
     - 0.007
     - 0.010
     - 0.014
   * - CF\ :sub:`2`
     - 0.002
     - 0.001
     - –
     - 0.001
   * - cyanogen
     - 0.013
     - 0.013
     - 0.009
     - 0.008
   * - glyoxal
     - 0.004
     - 0.004
     - 0.013
     - 0.016
   * - acrolein
     - 0.008
     - 0.008
     - 0.020
     - 0.007
   * - pyrazine
     - 0.002
     - 0.002
     - 0.002
     - 0.002

The evGW values agree to 0.001 Å for every state. G0W0 differs more
because the QP equations are solved in full here and linearized in the
Turbomole implementation; for CF\ :sub:`2` the G0W0 optimization opens
the F–C–F angle until S\ :sub:`1` nearly meets S\ :sub:`0`. Bonds to
hydrogen come out about 0.006 Å too long, as also reported by Holzer.
For acetylene in its :math:`^1\Sigma_u^-` state (def2-TZVPP) the QP
response moves the C–C minimum by +0.003 Å (G0W0) and +0.006 Å (evGW)
with the exact self-energy, in agreement with [Holzer:2026]_.


Embedding in frozen multipoles
******************************

With external multipoles (``dftgwbse`` ``mpsfile``, or sites given to the
DFT engine), the ground-state forces, the excited-state forces and the QP
gradients include the interaction with the sites at frozen moments:
point charges, dipoles (including induced dipoles of polar sites) and
quadrupoles, and the Gaussian induced charges and dipoles of
charge-response sites. The electrons enter with the ground-state density
or the relaxed density of the state (the sites are a one-electron
operator), the nuclei through the classical nucleus-site energy. The
forces on the sites are stored as well (orb datasets ``site_forces`` and
``excited_state_site_forces``, one row per site in the order of the mps
file), so that atoms and sites together obey Newton's third law.

The derivative integrals are three-centre Coulomb integrals between the
AO pair and Gaussian charges, dipoles and quadrupoles at the sites; point
multipoles are represented by narrow Gaussians whose width is extrapolated
to zero (agreement with finite differences of the point-multipole energy
to :math:`10^{-9}`). Frozen moments means that the induced moments of a
polarizable environment do not respond to the displacement; for an
environment converged self-consistently with the QM density their
response does not contribute to the ground-state forces (the polarization
energy is stationary), but it does for excited states.

Homogeneous field and Ewald background
**************************************

A homogeneous field (``externalfield``) adds
:math:`-\mathbf{F}\cdot(\mathbf{r}-\mathbf{R}_\mathrm{com})` to the
one-electron Hamiltonian and
:math:`\sum_A Z_A\mathbf{F}\cdot(\mathbf{R}_A-\mathbf{R}_\mathrm{com})` to
the energy of the nuclei, so that :math:`\mathrm{d}E/\mathrm{d}\mathbf{F}`
is the total dipole and a nucleus feels :math:`-Z_A\mathbf{F}` (the
nuclear term was added with the forces; DFT total energies in a field
change by it, excitation energies do not). The periodic Ewald background
of ``qmmm`` (``ewaldregion``) acts on the electrons through its potential
on the integration grid and on the nuclei through the potential at the
nuclei. For both, the electronic part of the force is the derivative of
the exact integral :math:`-\int (\partial\rho/\partial R_A)\,\phi`,
evaluated on the integration grid; finite differences of the energy, which
for the Ewald background is itself a grid quadrature, agree with it to the
accuracy of the grid (:math:`2\times 10^{-5}` relative on the ``fine``
grid). The gradient of the Ewald potential at the nuclei is taken by
central differences (:math:`h=10^{-4}` bohr) in the same evaluation. In a
``qmmm`` job (methane in the periodic methane background of the
integration tests, with a Thole region) finite differences of the total
energy agree with the forces to :math:`2\times 10^{-5}` (ground state)
and :math:`5\times 10^{-4}` (S1, limited by the near-degenerate state);
the Ewald part alone, from the difference of runs with and without the
background, to :math:`5\times 10^{-5}`. The
background is frozen, including the periodic images of the QM molecule
itself: the copy of the QM segment it leaves out is the lattice image of
the classical background copy nearest to the QM centroid, so small changes
of the QM geometry leave the background potential unchanged.

QM/MM forces
************

``qmmm`` with ``qmregion.compute_forces true`` computes, after the
inter-region loop has converged, the forces on the QM atoms of the tracked
state (ground state, charged e/h states, BSE singlets and triplets,
perturbative QP levels) with the environment at its converged moments,
and writes them to the job output (``forces``, Hartree/Bohr). They are
the derivatives of the converged QM/MM energy if the environment has been
converged against the density whose energy derivative it is: the
ground-state density for the ground state, and the relaxed density for
excited states (``relaxed_density true``; a warning is given otherwise).
The total energy is then stationary in the induced moments, so their
response does not enter the forces.
With ``environment_polarization ground`` and with ``environment_screening``
the environment is converged against the ground state, and its response
to the excitation's geometry dependence is part of the gradient (below).

Checks (methane S1, G0W0, 3-21G, finite differences of the converged
total energy): a self-consistent polarizable site in a unit test agrees
to :math:`6\times 10^{-6}`; a QM/MM run with one Thole region to
:math:`7\times 10^{-4}` and with static regions only to
:math:`1\times 10^{-4}` (the finite-difference accuracy); with two Thole
regions and a static region to :math:`1\times 10^{-4}` (ground state) and
:math:`8\times 10^{-4}` (S1).
With the relaxed density replaced by the unrelaxed one in the unit test
the forces are off by 7%.

Environment polarized by the ground state
*****************************************

With ``qmregion.environment_polarization ground`` the polar regions are
converged against the ground-state density, and the energy of an excited
state (or QP level) is the ground-state QM/MM energy plus the excitation
energy. The induced dipoles then follow the geometry through the ground
state, and the excitation energy depends on them, so frozen moments are
not enough. With the dipoles as variables of their own, the Lagrangian
gives two additional pieces:

* the Z-vector equation of the relaxed density contains the reaction field
  of the environment: a density change induces dipoles
  :math:`\delta\mu = -A^{-1} t`, :math:`t_i` its coupling to a dipole at
  site :math:`i`, :math:`A` the Thole operator of all polar regions
  together, whose potential acts back on the electrons (in the response,
  its source and the energy-weighted density);
* the dipoles :math:`\lambda = \delta\mu[P]` that the relaxed density
  :math:`P` would induce interact with the ground-state charges (density
  and nuclei); the derivative of that interaction is added.

The couplings are the exact AO integrals the QM Hamiltonian uses (an
RI-fitted field is off by several percent at sites in the tails of the
density). Checks: a polarizable site converged against the ground state
with S1 on top agrees with finite differences to
:math:`2\times 10^{-5}` (without the two pieces 4%); in QM/MM (methane
S1, two Thole regions) :math:`2\times 10^{-4}`, and along the direction
of the correction :math:`2\times 10^{-4}` against :math:`2\times 10^{-3}`
without it. With charge-response segments (methane model fitted with
``crfit``): a charge-response molecule converged against the ground state
:math:`3\times 10^{-5}` (5% without), QM/MM :math:`1.6\times 10^{-4}`,
along the correction :math:`1.4\times 10^{-4}` against
:math:`1.2\times 10^{-3}`. All polar regions with the same ``exp_damp``
(several polar regions: the default, which also couples them).

Screened GW-BSE
***************

In the screened scheme (``environment_screening``) the QM/MM loop
converges the polar regions against the ground state, and one GW-BSE run
on top screens the excitation by the environment: the three-centre
integrals are dressed with :math:`(1+R)^{1/2}` on the auxiliary index,
:math:`R` the reaction field of the polar sites in the metric of the fit
functions, and the QP energies contain the static reaction-field
self-energy :math:`\Sigma_\mathrm{reac}`. With ``qmregion.compute_forces``
the forces are those of :math:`E_\mathrm{GS}^\mathrm{QM/MM} + \Omega`:
the ground-state forces at the converged moments plus the gradient of
:math:`\Omega`, which includes the response of the environment through
the ground state (as with ``environment_polarization ground``, above) and
the dependence of the screening on the geometry.
``DFTEngine::setScreeningEnvironment`` makes the BSE gradient use the
same dressed integrals. The energy depends on the auxiliary representation
only through the effective interaction between the fit functions,
:math:`u = V^{-1} + V^{-1} B_\mathrm{env} V^{-1}`
(:math:`B_\mathrm{env} = -C^T \chi C`, :math:`C` the fields and
potentials of the auxiliary functions at the sites, :math:`\chi` the
fixed response of the environment), so its derivative follows from that
of the gas-phase code: the three-centre integrals, the Coulomb metric, and
:math:`B_\mathrm{env}`, whose sources move with the atoms
(``EnvironmentScreening::KernelGradient``; derivatives with respect to the
site positions by central differences). The static self-energies are not
functions of the dressed integrals: :math:`(1-a_x)\Sigma_x` is with the
bare :math:`V^{-1}`, :math:`\Sigma_\mathrm{reac}` with
:math:`V^{-1} B_\mathrm{env} V^{-1}`; the response of the QP corrections
(G0W0/evGW, ``ppm`` or ``exact``) differentiates both kernels explicitly.
Without ``include_kreac`` the BSE exchange is with the bare
:math:`V^{-1}` as well, and is differentiated the same way; the direct
term is unchanged, since :math:`W_\mathrm{tot} = (u^{-1} + \epsilon - 1)^{-1}`
with the bare integrals gives the same interaction as
:math:`\epsilon^{-1}` with the dressed ones.

Checks: S1 of water in a fixed Thole environment (explicit and shell
sites, smeared) and in a charge-response molecule, QP corrections fixed:
finite differences :math:`3\times 10^{-6}` and :math:`3\times 10^{-5}`
(the screening moves the projected gradient by 40%). ``qmmm`` (methane
S1, G0W0, two Thole regions, one of them a shell with
:math:`\epsilon = 4`): :math:`7\times 10^{-5}` along the forces,
:math:`8\times 10^{-4}` along a random direction (absolute
:math:`4\times 10^{-8}` Ha); with the QP corrections kept fixed the
forces are off by 7% and 28% along those directions, without the
response of :math:`\Sigma_\mathrm{reac}` by 2% and 18%. Without
``include_kreac``: water in the Thole environment
:math:`2\times 10^{-5}`; ``qmmm`` along the forces
:math:`1\times 10^{-5}` (Richardson extrapolation of the steps
:math:`5\times 10^{-4}` and :math:`10^{-3}` Å, which differ from the
forces by :math:`2.4\times 10^{-4}` and :math:`9.9\times 10^{-4}`, the
quadratic truncation of the near-degenerate S1 of methane),
:math:`6\times 10^{-4}` along a random direction.

Geometry optimization in QM/MM
******************************

``qmmm`` itself does not move atoms; an optimizer drives it by running one
job per geometry. With ``compute_forces`` the job output also holds the QM
coordinates the forces belong to, in the frame of the job (``geometry``,
Angstrom), and the job directory holds the final orbitals
(``qmregion_final.orb``). The next step gives the new coordinates as an xyz
file in that frame (``qmregion.geometry``), which replaces the segment
mapping: the mapping fits each fragment onto its MD frame and would undo
rigid-body motion and change the geometry between fragments. The orbitals of
the previous step can start the SCF (``qmregion.guess``; they are
orthogonalized in the new geometry), and the induced moments of the
previous step the polar regions (``polarregion.guess``, a checkpoint of
that job; segments are matched by id, and a checkpoint that does not fit
the region is ignored with a warning).

The environment is frozen at its topology: the region partition follows
the segment positions of the topology, and the Ewald background, including
the periodic images of the QM molecule, stays as converged (the copy of the
QM segment it leaves out is fixed by the topology). The polar regions are
re-converged in every job, so with ``relaxed_density`` the energy is
stationary in their induced moments and the forces are its derivatives.
If the geometry changes much, the multipoles of the QM segment can be
refitted at the optimized geometry and the background recomputed (usually
one such cycle).

Check (the S1 setup above, H1 displaced by :math:`\pm 5\times 10^{-4}` Å
through ``qmregion.geometry``, guess from the reference job): finite
differences agree with the forces to :math:`8\times 10^{-4}`; the
orbital guess halves the number of iterations of the first SCF, and the
induced moments of the reference job reduce the QM/MM iterations from 19
to 14 at unchanged energy and forces.

``xtp_qmmm_optimize`` (with ``pyxtp``) runs this loop with an ASE optimizer
(BFGS, FIRE or Sella)::

    xtp_qmmm_optimize -o qmmm.xml -f state.hdf5 -t 4 --fmax 0.05

With ``--optimizer Sella`` the QM region is optimized in internal
coordinates plus explicit translations and rotations of the whole region:
in a fixed environment the energy is not invariant under rigid-body motion,
which Sella projects out by default (as for a free molecule) and would then
stop with a net force and torque left. Sella converges the soft modes of a
solute in a cage (torsions, out-of-plane modes, rigid-body motion) much
faster than BFGS in Cartesian coordinates.
A frozen cage from a classical MD snapshot need not fit the QM model (other
charges, polarization, the excited state), and the QM region can then
drift and rotate in it for many steps (anisole S1 in water: 0.28 Å and
13° within 17 steps, 1 meV per step). ``--freeze-rigid-body`` (Sella)
keeps the position and orientation of the QM region and relaxes only its
internal coordinates; convergence then refers to the forces without net
force and torque.

It takes the qmmm options and job file of an existing setup (the first job,
or ``--job-id``), passes the orbitals and induced moments of each step on
to the next, runs every step in its own directory under ``--workdir``
(the entries of the setup directory are linked there; the state file is
only read), and writes the trajectory and the optimized geometry in the job
frame (``optimized.xyz``, for ``qmregion.geometry`` of later jobs). Atoms
can be kept in place with ``--fix``. In Python,
``pyxtp.qmmm_optimize.QMMMCalculator`` is the ASE calculator behind it.
A job that fails (e.g. a QM/MM loop that does not converge in
``max_iterations`` once the state changes character) stops the
optimization with the step directory to look at; check that the tracked
state stays the intended one (the summary prints it per step).

Short-range terms: the qmmm energy has no repulsion between QM and
environment atoms (the polar and static regions are point multipoles,
point or Gaussian induced moments), and the xtp DFT has no dispersion.
Single points are fine with that; an optimization is not: in water, the
hydroxyl H of phenol moves onto the O of a neighbouring water and the O-H
bond stretches (0.96 to 1.00 Å in a few steps, energy still falling). As
in standard QM/MM, ``--lj params.json`` (``pyxtp.qmmm_optimize.
LennardJones``) adds Lennard-Jones terms between QM and environment atoms,
typically with the parameters of the force field of the MD (JSON:
``"qm"`` element table or one ``[sigma, epsilon]`` per QM atom,
``"environment"`` element table, optional ``"cutoff"``; σ in Å, ε in
kJ/mol, Lorentz-Berthelot mixing; environment elements without an entry get
none). ``--d3 pbe0`` (``pyxtp.qmmm_optimize.D3Dispersion``, needs the python
package ``dftd3`` of simple-dftd3) adds Grimme's D3 with the damping
parameters of the given functional (``--d3-damping``, default
Becke-Johnson; three-body term with ``--d3-atm``) within the QM region and,
without ``--lj``, also with the environment atoms within ``--d3-cutoff``
(15 Å) of a QM atom (``--d3-qm-only`` for the QM region alone; with
``--lj`` the dispersion between QM and environment is that of LJ). The
environment atoms are the sites with an element name in the checkpoint of
the first job; they stay where they are, each term keeps the atoms within
its cutoff at the start, and the D3 energy with the environment is
:math:`E_\mathrm{D3}(\mathrm{QM+env}) - E_\mathrm{D3}(\mathrm{env})`, so
the forces on the QM atoms are its exact derivatives. These terms enter
the geometry only: the electronic structure (and the excitation energy) is
that of the job, and the same terms are added for ground and excited
states.

Limitations
***********

* Closed-shell ground states; singlet and triplet BSE states.
* Only the xtp DFT package, with an auxiliary basis (RI).
* Embedding in frozen external multipoles, homogeneous fields and the
  frozen Ewald background, and polar regions converged against the state
  (``relaxed_density``) or the ground state (``environment_polarization
  ground``; see above), and the screened scheme
  (``environment_screening``); no external
  densities or potential matrices, no truncated bases.
* No off-diagonal QP Hamiltonian in the BSE (``use_Hqp_offdiag``), no
  QSGW, no contour deformation (QP corrections are then kept fixed), no
  dynamical screening of the BSE kernel.
* No effective core potentials.
* Degenerate KS levels are assumed to be symmetry partners.
