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
with equal numbers of :math:`\alpha` and :math:`\beta` electrons, and an
unrestricted-Kohn-Sham (UKS) class for everything else. The two share the same set of
convergence options and the same second-order (TRAH) solver described below.

Each SCF iteration builds a Fock matrix guess by one of three mutually-exclusive
methods, in order of preference:

- **A-DIIS** (energy-based DIIS) [Hu:2010]_, used while the SCF error is still
  relatively large, and switching to
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
mechanisms in ORCA's own SCF converger.

Second-order (TRAH) steps
=========================

DIIS-family extrapolation is well known to struggle for open-shell systems with a
small HOMO-LUMO gap, or for systems whose physically correct electronic configuration
is not the solution DIIS happens to converge to; DIIS can also come to rest on a
stationary point of the energy that is not a minimum. For these cases the engine can
minimize the energy directly in orbital-rotation space with trust-region
augmented-Hessian (TRAH) steps [Helmich-Paris:2021]_. The option
``xtpdft.convergence.second_order`` selects when:

- ``auto`` (default): DIIS-family extrapolation and mixing as above, switching to TRAH
  steps once (A)DIIS has failed five times in a row, or once, after at least 50
  iterations, the DIIS error over the last 10 iterations has on average shrunk by less
  than a factor 1.125 per iteration. The second criterion is checked in every
  iteration, also when the error never falls below ``ADIIS_start`` and (A)DIIS is never
  tried.
- ``always``: TRAH steps from the first iteration.
- ``off``: never.

Once switched on, TRAH steps are used for the rest of the SCF. They are taken for
closed-shell RKS and for UKS calculations; restricted open-shell and fractional
occupations are not covered.

**Parametrization.** The orbitals are rotated as :math:`\mathbf{C}' = \mathbf{C}\,
\exp(\boldsymbol{\kappa})` with :math:`\kappa_{ia} = x_{ia}` and :math:`\kappa_{ai} =
-x_{ia}` (:math:`i` occupied, :math:`a` virtual; one block per spin channel for UKS).
The exponential is evaluated exactly from the singular value decomposition of the
occupied-virtual block, so orthonormal orbitals stay orthonormal. With :math:`F` the
Fock matrix in the current orbitals, the energy gradient is

.. math::

    \frac{\partial E}{\partial x_{ia}} = -4 F_{ia} \quad \textrm{(RKS)}, \qquad
    \frac{\partial E}{\partial x^\sigma_{ia}} = -2 F^\sigma_{ia} \quad \textrm{(UKS)}.

**Hessian.** The Hessian is never stored; the solver only needs its product with trial
vectors, which is computed analytically. For RKS, with :math:`X_{ai} = x_{ia}`,

.. math::

    (\mathbf{H}\mathbf{x}) = 4 \left[ F_{vv} X - X F_{oo}
      + C_v^T\, G[\delta P]\, C_o \right], \qquad
    \delta P = 2 \left( C_v X C_o^T + C_o X^T C_v^T \right),

where :math:`G[\delta P] = J[\delta P] + \tfrac{1}{2} a_x K[\delta P] + \delta V_{xc}`
is the linear response of the Coulomb, exact-exchange and exchange-correlation
potentials, the latter from the second functional derivative :math:`f_{xc}` of the
exchange-correlation energy (the same response kernel as in the static
coupled-perturbed Kohn-Sham calculation of the polarizability). For UKS,

.. math::

    (\mathbf{H}\mathbf{x})_\sigma = 2 \left[ F^\sigma_{vv} X_\sigma - X_\sigma F^\sigma_{oo}
      + C^{\sigma\,T}_v\, G_\sigma\, C^\sigma_o \right], \qquad
    G_\sigma = J[\delta P_\alpha + \delta P_\beta] + a_x K[\delta P_\sigma]
      + \delta V_{xc,\sigma}[\delta P_\alpha, \delta P_\beta],

with :math:`\delta P_\sigma = C^\sigma_v X_\sigma C^{\sigma\,T}_o + \textrm{h.c.}` and the
spin-polarized :math:`f_{xc}`, so that the two spin channels are coupled through the
Coulomb and exchange-correlation response. Both work with RI and with four-center
integrals, for LDA, GGA and global hybrid functionals. If a functional has no second
derivatives in libxc, the UKS Hessian products fall back to central finite differences
of the Fock matrix (two Fock builds per product); the CDFT constraint potential is part
of :math:`F` and has no response term.

**Step.** Each step minimizes the quadratic model :math:`\mathbf{g}^T\mathbf{x} +
\tfrac{1}{2}\mathbf{x}^T\mathbf{H}\mathbf{x}` within a trust radius :math:`R`. The
lowest eigenpair :math:`(\mu, [v_0; \mathbf{v}])` of the scaled augmented Hessian

.. math::

    \begin{pmatrix} 0 & \alpha\, \mathbf{g}^T \\ \alpha\, \mathbf{g} & \mathbf{H} \end{pmatrix}

gives the step :math:`\mathbf{x} = \mathbf{v} / (\alpha v_0)`, which satisfies
:math:`(\mathbf{H} - \mu)\mathbf{x} = -\mathbf{g}`. It is found with a Davidson
iteration (preconditioned by the diagonal :math:`4(F_{aa} - F_{ii})`, respectively
:math:`2(F^\sigma_{aa} - F^\sigma_{ii})`). Since only the first row and column of the
projected matrix depend on :math:`\alpha`, all values of :math:`\alpha` are handled in
the same subspace without further Hessian products: :math:`\alpha = 1` is used if its
step lies inside the trust radius, otherwise :math:`\alpha` is increased by bisection in
:math:`\log\alpha` until :math:`|\mathbf{x}|` lies between :math:`0.8R` and :math:`R`.
The micro-iterations stop when the residual :math:`|(\mathbf{H}-\mu)\mathbf{x} +
\mathbf{g}|` falls below :math:`\min(0.1, |\mathbf{g}|)\,|\mathbf{g}|`, which makes the
SCF converge quadratically close to the solution; ``davidson_max_iter`` limits their
number.

**Trust radius.** A step is checked in the next iteration, once the energy of its
density is known, by the ratio :math:`r` of actual to predicted energy change
[Helmich-Paris:2021]_. It is rejected if :math:`r < 0`; otherwise the radius is reduced
by 0.7 for :math:`r \le 0.25` and increased by 1.2 for :math:`r > 0.75` (initial radius
0.5, at most 1.0). After a rejection, the energy along the step,
:math:`E(t) = E_0 + (\mathbf{g}^T\mathbf{x})\, t + c\, t^2`, is fitted to the actual
change, and the new radius is set to its minimum :math:`t^*|\mathbf{x}|`, with
:math:`t^*` kept between 0.1 and 0.7. The rejected step is undone and repeated at once
from the stored Fock matrix and Davidson subspace, so it costs neither an additional
Fock build nor, usually, any Hessian products. If the radius falls below
:math:`3\times10^{-3}` without an accepted step, the SCF continues without TRAH steps.
If the density handed to the converger is not the density of the current orbitals
(after a damped iteration, or a guess density), the orbital density is returned first
and the step is taken in the following iteration.

**What to expect.** TRAH steps converge to a stationary point by lowering the energy in
every accepted step. In difficult cases (for example a stretched O\ :sub:`2` triplet or
an iron(II) complex) this is often a lower state than the one DIIS converges to, which
can be a saddle point. The end point is not guaranteed to be a minimum, however:
directions of negative curvature that break a symmetry of the starting orbitals carry
no gradient and are not followed, since no stability analysis is performed. Per
iteration a TRAH step costs more than a DIIS step (several Hessian products), so
``auto`` remains the default; ``always`` is useful for systems that are known to be
difficult.

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
DIIS-family, direct-minimization, or CDFT) produced that reference. The methods,
numerics and settings of the GW and BSE steps -- level windows, self-energy models,
the solution of the quasiparticle equation, evGW self-consistency and its pitfalls --
are described in detail in :doc:`gwbse_numerics`.
