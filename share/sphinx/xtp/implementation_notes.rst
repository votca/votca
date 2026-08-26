Implementation Notes
#####################

This page collects practical, code-level notes on the internal DFT, CDFT, and GW-BSE
implementations that go beyond the level of detail in :doc:`elstruc_theory`. It is aimed
at developers and users who need to understand *why* a calculation behaves the way it
does, not just the underlying theory.

DFT: self-consistent field convergence
***************************************

The internal DFT engine dispatches to one of two convergence-acceleration classes
depending on the system: a restricted-Kohn-Sham (RKS) class for closed-shell systems
with equal numbers of :math:`\alpha` and :math:`\beta` electrons, and a separate,
unrestricted-Kohn-Sham (UKS) class for everything else. The two share the same set of
convergence options, but are otherwise independent implementations; nothing described
below for the UKS path affects RKS calculations.

Each SCF iteration builds a Fock matrix guess by one of three mutually-exclusive
methods, in order of preference:

- **A-DIIS** (energy-based DIIS), used while the SCF error is still relatively large,
  and switching to
- **plain (commutator) DIIS** [Pulay:1980]_ once the error is small, both extrapolating
  a new Fock matrix guess from a short history of previous Fock and density matrices,
  or
- **linear mixing** of the current density with that of the previous iteration, used
  whenever DIIS-family extrapolation is judged to be unreliable for the current
  iteration.

Independently of which of these is chosen, two further, always-available corrections
are applied to the resulting Fock matrix before it is diagonalized: level shifting of
the virtual-orbital block, and (for the mixing branch specifically) an *adaptive*
damping parameter that ramps from a base value toward a configurable ceiling as
consecutive convergence failures accumulate, and relaxes back toward the base value as
readily as it climbed once the SCF starts succeeding again. Both the damping ramp and
a trailing-average stall detector (see below) are modeled directly on the equivalent
mechanisms in ORCA's own SCF converger, whose behavior was compared against directly
during development of this feature.

Direct-minimization fallback
=============================

DIIS-family extrapolation is well known to struggle for open-shell systems with a
small HOMO-LUMO gap, or for systems whose physically correct electronic configuration
is not obviously the lowest-energy solution DIIS would otherwise converge to. When
either

- DIIS-family extrapolation has failed for several consecutive iterations, or
- the SCF error is failing to improve by more than a set tolerance, on average, over a
  trailing window of iterations (evaluated only after a minimum number of total
  iterations, so that this check is not triggered spuriously early in a calculation),

the UKS engine falls back to direct energy minimization in orbital-rotation space,
rather than continuing to rely on DIIS-family extrapolation or plain mixing alone. This
mirrors the automatic switch to a Trust-Region Augmented-Hessian (TRAH) converger found
in mature quantum-chemistry codes such as ORCA.

The fallback solves the augmented-Hessian eigenvalue problem of [Helmich-Paris:2021]_,

.. math::

    \begin{pmatrix} 0 & \alpha\, \mathbf{g}^T \\ \alpha\, \mathbf{g} & \mathbf{H} \end{pmatrix}
    \begin{pmatrix} 1 \\ \boldsymbol{\kappa}(\alpha) \end{pmatrix}
    = \mu \begin{pmatrix} 1 \\ \boldsymbol{\kappa}(\alpha) \end{pmatrix} ,

for its lowest eigenpair, which simultaneously yields a level shift :math:`\mu` and an
orbital-rotation vector :math:`\boldsymbol{\kappa}`. Unlike a fixed, diagonal Hessian
approximation, this correctly handles the indefinite Hessian directions that occur near
genuine electronic near-degeneracies, at the cost of needing an iterative solve (via a
Davidson eigensolver) rather than a single diagonalization. The scalar :math:`\alpha` is
adjusted by bisection so that the resulting step stays within an adaptive trust radius,
following the same step-control logic (accept, shrink, or grow the trust radius
depending on how well the step's predicted energy change matches the energy change
actually observed) described in [Helmich-Paris:2021]_.

The Hessian-vector product required by the Davidson solver ("sigma vector") is
evaluated by a *central* finite difference of the orbital gradient, rather than the
exact, analytic coupled-perturbed response. This is a deliberate simplification:
implementing the exact response would require a second, coupled-perturbed Fock build
using the second functional derivative of the exchange-correlation energy, which has
not been implemented. A central (rather than one-sided) finite difference is used
specifically because a one-sided difference does not preserve the exact symmetry of the
true Hessian, which was found in practice to break the Davidson solver's own
convergence.

For open-shell systems, the :math:`\alpha` and :math:`\beta` spin channels can either be
treated independently -- solving two separate augmented-Hessian problems, one per
channel, holding the other channel's density fixed -- or as a single, coupled problem
over the combined rotation space of both channels together. The coupled treatment is
generally more robust, since it accounts for the coupling between the two channels
through the shared Coulomb potential and the exchange-correlation kernel's cross-spin
terms; the independent treatment is cheaper per step, since each of the two smaller
eigenvalue problems is faster to solve than the one, larger combined problem. Both are
available; the coupled treatment is preferred whenever the necessary Fock-matrix-build
callback has been supplied.

Constrained DFT
***************

CDFT is implemented as two nested loops around the ordinary UKS SCF machinery described
above: an inner loop that runs an ordinary UKS SCF for a fixed value of the constraint
potential :math:`V_c`, and an outer bisection loop that adjusts :math:`V_c` until the
constrained population, evaluated once the inner SCF has converged, matches the
requested target to within a set tolerance. The constraint potential itself,
:math:`V_c\, w_c(\mathbf{r})`, is added directly to the Fock matrix inside the inner
SCF, so it participates in DIIS-family extrapolation exactly like the rest of the Fock
matrix; no special handling is needed there.

Each new inner SCF is warm-started from the previous one's converged orbitals, since
neighboring values of :math:`V_c` along the bisection typically have similar,
neighboring solutions. This is generally beneficial, but is not guaranteed to be: if a
bisection step is large (for example, during the initial bracket-expansion phase, before
the root has been localized to a narrow interval) and the electronic structure changes
character across that interval, warm-starting from the previous, qualitatively
different solution can bias the new inner SCF toward the wrong electronic state rather
than simply providing a mildly imperfect starting point.

Hirshfeld weight function and forces
=====================================

The Hirshfeld partitioning of :eq:`equ:hirshfeld_weight` is evaluated on the same
real-space integration grid used for the exchange-correlation energy, using free-atom
reference densities computed once per element at the start of the calculation. Because
this weight function depends on the nuclear positions, the CDFT contribution to the
nuclear forces requires several terms beyond the ordinary Hellmann-Feynman force,
following [Oberhofer:2009]_:

- the derivative of the weight function itself with respect to the nuclear positions,
- the derivative of the grid point weights (from the space-partitioning scheme
  underlying the numerical integration) with respect to the nuclear positions,
- a Pulay-type term from the basis functions themselves moving with the nuclei, and
- a term from the translation of the integration grid points as the nuclei move.

All four terms are evaluated explicitly and summed, and have been validated against
finite-difference forces on constrained systems.

GW-BSE
******

The :math:`GW` and Bethe-Salpeter equation implementations, described theoretically in
:doc:`elstruc_theory`, are built on top of the ground-state DFT orbitals and Fock matrix
produced by the machinery above -- either from the internal DFT engine or from the
ORCA interface -- and do not themselves depend on which SCF convergence path (ordinary
DIIS-family, direct-minimization, or CDFT) produced that reference. Consult the code
reference for implementation-level detail on the GW self-energy and BSE Hamiltonian
construction; this page focuses specifically on the ground-state SCF machinery that
CDFT and difficult open-shell systems most directly exercise.
