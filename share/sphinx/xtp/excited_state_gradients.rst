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
* The restrictions of the excited-state forces apply (closed shell, RI, no
  embedding, G0W0 or evGW with ``ppm`` or ``exact``; otherwise the
  gradient of the KS energy with fixed corrections is returned).

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


Limitations
***********

* Closed-shell ground states; singlet and triplet BSE states.
* Only the xtp DFT package, with an auxiliary basis (RI).
* No embedding or environment screening (QM/MM), no truncated bases.
* No off-diagonal QP Hamiltonian in the BSE (``use_Hqp_offdiag``), no
  QSGW, no contour deformation (QP corrections are then kept fixed), no
  dynamical screening of the BSE kernel.
* No effective core potentials.
* Degenerate KS levels are assumed to be symmetry partners.
