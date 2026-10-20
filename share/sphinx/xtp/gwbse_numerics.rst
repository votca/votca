.. _gwbse_numerics:

GW-BSE in Practice: Methods, Numerics and Settings
##################################################

:doc:`elstruc_theory` introduces the quasiparticle (QP) equation and the
Bethe-Salpeter equation (BSE) [Hedin:1965]_ [Hybertsen:1986]_
[Rohlfing:2000]_. This chapter describes how VOTCA-XTP solves them: which
levels enter where, how the self-energy is evaluated, how the QP equation is
solved for each level, how the eigenvalue self-consistency
(evGW) is iterated, and how the BSE is set up and solved. It also explains
the numerical difficulties one runs into in practice -- several solutions of
the QP equation, levels that hop between them, self-consistency loops that
stall -- and what the settings do about them. All options mentioned below are
in the ``gwbse`` section of the options file; their defaults are those of the
current version.

.. contents:: Contents of this chapter
   :local:
   :depth: 2


Overview of a GW-BSE calculation
********************************

A ``dftgwbse`` or QM/MM GW-BSE job runs through the following steps:

#. **Level windows.** From the DFT orbitals, three ranges of molecular
   orbitals (MOs) are chosen: the levels that enter the screening (RPA
   window), the levels for which QP energies are computed (GW window) and
   the levels from which the electron-hole pairs of the BSE are built (BSE
   window). See :ref:`gwbse_windows`.
#. **Three-centre integrals.** With the resolution of identity (RI) in the
   auxiliary basis, the Coulomb interaction between orbital products is
   expressed by three-centre integrals :math:`M^P_{mn}`, with
   :math:`m` in the GW (or BSE) window, :math:`n` in the RPA window and
   :math:`P` an auxiliary function, in an orthonormalised auxiliary basis in
   which the bare Coulomb interaction is :math:`v_{mn,kl}=\sum_P
   M^P_{mn}M^P_{kl}`.
#. **Exchange self-energy** :math:`\Sigma_\mathrm{x}`, computed once.
#. **Screening.** The screened interaction :math:`W` from the
   random-phase approximation (RPA), either through a plasmon-pole model
   (PPM) or exactly through the RPA excitation spectrum. See
   :ref:`gwbse_sigma`.
#. **QP equation** for every level of the GW window, see :ref:`gwbse_qp`.
#. **Self-consistency** (evGW): the QP energies replace the DFT energies in
   :math:`G` and :math:`W`, and steps 4 and 5 are repeated until the energies
   no longer change, see :ref:`gwbse_evgw`.
#. **QP Hamiltonian.** The off-diagonal elements of :math:`\Sigma` give
   the QP Hamiltonian :math:`H^\mathrm{QP}` and, by diagonalising it, the
   "diagonalized" QP energies (DQP) next to the perturbative ones (PQP).
#. **BSE** with static screening for singlets and/or triplets, optionally
   followed by a perturbative correction for dynamical screening, see
   :ref:`gwbse_bse`.

The log file reports each step with time stamps. At the default log level
(``--verbose1`` for more, ``-v`` for debug output) it contains the level
ranges, a summary of every screening build (number of active PPM modes, or
RPA correlation energy and lowest RPA poles for exact screening), the evGW
iteration summary, the QP and BSE energies and a timing
summary of the GW part.


.. _gwbse_windows:

Level windows
*************

Ranges by index
===============

The option ``ranges`` sets the windows by orbital index:

``default``
    all levels in the RPA; GW and BSE from the lowest level up to
    :math:`3\,\mathrm{HOMO}+1`, i.e. twice as many virtual as occupied levels;
``factor``
    ``rpamax``, ``qpmin``, ``qpmax``, ``bsemin``, ``bsemax`` as fractions of
    the number of levels or of occupied levels;
``explicit``
    the same five values as orbital indices;
``full``
    all levels everywhere.

The ``default`` choice is generous. For larger molecules its upper end lies
20 eV or more above the LUMO (C60: level 538, about 24 eV above the LUMO;
a D5M dimer: level 994). These levels are expensive -- every QP iteration
solves the QP equation for each of them -- and they lie in the energy range
where the QP picture becomes ambiguous (see :ref:`gwbse_satellites`).

Energy caps
===========

The options ``qp_emax``, ``qp_emin``, ``bse_emax`` and ``bse_emin`` (in eV;
0 means off, the default) shrink the windows set by ``ranges`` to the levels
within the given energy above the DFT LUMO (``*_emax``) or below the DFT HOMO
(``*_emin``)::

    <qp_emax>15</qp_emax>
    <bse_emax>15</bse_emax>

The caps only ever remove levels. A set of levels within 5·10\ :sup:`-4` Ha
of the edge level is kept whole, so that a degenerate set is never cut. For
open-shell (UKS) references a level is kept if it is within the cap in
either spin channel. The log reports every cap, for example
``GW window top capped by qp_emax = 20 eV: 13 -> 11 (2 levels removed)``.

Which cap is appropriate is a physical choice. The low-lying excitations of
organic molecules are usually built from transitions into virtual levels
within 10–15 eV of the LUMO; the trouble with satellites described below
typically starts 20 eV or more above it.

Core levels
===========

``ignore_corelevels`` removes the chemical core levels (counted per element)
from parts of the calculation:

``RPA``
    from the screening, the GW window and the BSE (frozen core);
``GW``
    from the GW window and the BSE, but they still screen;
``BSE``
    from the BSE only;
``none``
    (default) all levels everywhere.

With ``RPA`` the three-centre integrals start at the first valence level,
so the exchange self-energy would only sum over the occupied valence levels,
while :math:`V_\mathrm{xc}` contains the core density. ``core_exchange``
(default ``true``) adds the missing exchange with the core levels,

.. math::

   \Delta\Sigma^\mathrm{x}_{mn} = -\sum_{c\in\mathrm{core}} (mc|cn),

computed once from three-centre integrals filled for this purpose. With it,
frozen-core GW reproduces all-electron GW closely: for a QM9 molecule
(def2-TZVP, 9 core levels, evGW with exact :math:`\Sigma`) the exchange
self-energies agree to all printed digits, the QP energies within 0.01 eV
and the lowest singlet within 3 meV, at 30 % lower cost. Without the core
exchange the occupied QP energies are too high by 0.6 eV on average (up to
1.4 eV) and the virtual ones by 0.25 eV. ``core_exchange=false`` restores
that behaviour for comparison with old results.

Levels outside the GW window
============================

Levels in the RPA (or BSE) window but outside the GW window have no QP
energy of their own. In evGW their energies in :math:`G` and :math:`W` are
the DFT energies shifted according to ``out_of_window_shift``:

``max`` (default)
    by the largest QP correction of the occupied (for levels below the
    window) or virtual (above) window levels;
``nearest``
    by the correction of the window level next to the boundary;
``linear``
    by a linear fit of the corrections of the window levels next to the
    boundary (a fraction ``out_of_window_fit_fraction`` of them, at least
    two), extrapolated at most over the energy span of the fit.

The shift is updated in every evGW iteration. The BSE uses the same shifted
energies, without off-diagonal QP elements, for virtual levels above the GW
window (if ``bse_emax`` or ``bsemax`` reaches further than ``qp_emax`` or
``qpmax``) and occupied levels below it. The further the BSE extends beyond
the GW window, the more the choice of shift matters; ``linear`` follows the
trend of the corrections near the edge.


.. _gwbse_sigma:

The correlation self-energy
***************************

Pole structure
==============

All self-energy models used for the QP equation give the diagonal
correlation self-energy of a level :math:`n` as a sum over poles,

.. math::
   :label: equ:gwbse:polesum

   \Sigma^\mathrm{c}_{nn}(\omega) = \sum_{m}\sum_{s}
     \frac{a_{nms}\,(\omega - z_{ms})}{(\omega - z_{ms})^2 + \eta^2},
   \qquad
   z_{ms} = \begin{cases} \varepsilon_m - \Omega_s & m\ \text{occupied}\\
                          \varepsilon_m + \Omega_s & m\ \text{virtual}
            \end{cases}

with :math:`m` running over the RPA window, :math:`s` over the modes of the
screening with frequencies :math:`\Omega_s`, non-negative weights
:math:`a_{nms}` and the broadening ``eta`` (default :math:`10^{-3}` Ha). The
poles :math:`z_{ms}` combine an orbital energy in :math:`G` with an
excitation energy of the screening. This structure determines how the QP
equation behaves (:ref:`gwbse_qp`).

Plasmon-pole model (``sigma_integrator=ppm``, default)
======================================================

The PPM represents the inverse dielectric matrix in the auxiliary basis by
one pole per eigenmode, fitted at a real and an imaginary frequency in the
spirit of [Godby:1989]_. The dielectric matrix is computed in the RPA at two
frequencies, :math:`\omega=0` and the imaginary frequency :math:`\omega =
0.5\,i` Ha. Its eigenvectors at :math:`\omega=0` define the modes, their
eigenvalues the weights, and the second frequency fixes the pole frequency
:math:`\Omega_s` of each mode. Modes with weight below :math:`10^{-5}` are
dropped; modes whose fit gives a negative squared frequency are counted and
reported in the log (their absolute value is used). The three-centre
integrals are then rotated into the frame of the modes, so that
:math:`a_{nms}\propto w_s\,\Omega_s\,|M^s_{nm}|^2`.

Cost per screening build: two dielectric matrices (:math:`N_\mathrm{aux}^2
N_\mathrm{occ}N_\mathrm{virt}`), one eigendecomposition
(:math:`N_\mathrm{aux}^3`) and the rotation of the integrals
(:math:`N_\mathrm{GW}N_\mathrm{RPA}N_\mathrm{aux}^2`). For large systems
the rotation dominates: for the D5M dimer (3195 orbitals, 7949 auxiliary
functions, 995 GW levels) it takes about 600 s per build on 28 cores, of a
total of about 900 s.

Exact RPA screening (``sigma_integrator=exact``)
================================================

The poles of :math:`W` are the RPA excitation energies, obtained from the
symmetric pair eigenproblem

.. math::

   (A-B)^{1/2}(A+B)(A-B)^{1/2}\,\mathbf{X}_s = \Omega_s^2\,\mathbf{X}_s ,

of dimension :math:`N_\mathrm{occ}N_\mathrm{virt}` (the number of
occupied-virtual pairs in the RPA window), with :math:`A-B` diagonal for the
RPA. The residues follow from the eigenvectors through the screening modes
(:math:`N^5` rather than :math:`N^6`), and :eq:`equ:gwbse:polesum` holds with
:math:`s` running over all RPA excitations. The eigendecomposition scales as
:math:`(N_\mathrm{occ}N_\mathrm{virt})^3` and dominates for anything but
small molecules (QM9 molecule with def2-TZVP: 7775 pairs, about 17 s per
build); it limits the method to molecules of a few dozen atoms. In exchange
there is no model error.

Contour deformation (``sigma_integrator=cda``)
==============================================

The frequency integral of :math:`G W` is deformed into residues of the
poles of :math:`G` enclosed by the contour plus an integral along the
imaginary axis, evaluated with a Gauss quadrature (``quadrature_scheme``,
``quadrature_order``); ``alpha`` smooths the residue term. The imaginary-axis
and tail terms are computed once per screening build. CDA has no
off-diagonal elements of :math:`\Sigma`.

Batched evaluation and far-pole interpolation
=============================================

The QP search (next section) evaluates :math:`\Sigma^\mathrm{c}_{nn}` of one
level at many frequencies at once -- the 60 or so points of the coarse scan
of a window of :math:`\pm0.75` Ha. For PPM and exact screening these batches
use the pole structure of :eq:`equ:gwbse:polesum`: poles far from the
frequency range :math:`[c-h, c+h]` are summed at Chebyshev nodes of the range
and interpolated, poles close to it are summed at every frequency. For a pole
at distance :math:`d=u\,h` from the centre, the kernel is analytic inside the
Bernstein ellipse with :math:`\rho = u+\sqrt{u^2-1}`, and :math:`n`-node
interpolation converges like :math:`\rho^{-n}` [Trefethen:2013]_. The poles are grouped by
:math:`u` (48 nodes for :math:`2\le u<3`, 32 for :math:`3\le u<5`, 24 for
:math:`5\le u<10`, 16 beyond), so that the interpolation error stays below
:math:`10^{-20}` of each pole's contribution; poles with :math:`u<2` are
summed directly. The result agrees with the term-by-term sum to rounding
(:math:`10^{-13}` relative in the unit tests). The split is used from 24
frequencies on. For a dimer-sized system it makes a 63-point batch 2.3–2.7
times faster, since only 4–10 % of the poles are near the window.

Single frequencies (the refinement steps of the root search) and
derivatives are summed term by term. Note that near a pole, where
:math:`|\omega-z|` is of the order of :math:`\eta`, rounding of
:math:`\omega - z` alone changes a term by up to :math:`a\,\epsilon\,|\omega|
/\eta^2`; values at the same frequency computed in different ways can
therefore differ at the :math:`10^{-12}`–:math:`10^{-11}` relative level,
depending on compiler and vectorisation.

Off-diagonal elements
=====================

The off-diagonal elements needed for :math:`H^\mathrm{QP}` are
evaluated as blocked matrix products for all pairs at once (PPM, exact).
They are only needed for the diagonalised QP energies and for the BSE with
``use_Hqp_offdiag``; for 995 levels (about :math:`5\cdot10^5` pairs) they take
of the order of ten minutes with PPM.


.. _gwbse_qp:

Solving the quasiparticle equation
**********************************

The QP equation and its roots
=============================

For each level :math:`n` of the GW window the perturbative QP energy is a
root of

.. math::
   :label: equ:gwbse:qpfunc

   f_n(\omega) = \varepsilon_n^\mathrm{KS} + \Sigma^\mathrm{x}_{nn}
                 - V^\mathrm{xc}_{nn} + \Sigma^\mathrm{c}_{nn}(\omega) - \omega = 0 ,

and its quasiparticle weight (spectral weight, renormalisation factor) is

.. math::

   Z_n = \left(1 - \left.\frac{\partial\Sigma^\mathrm{c}_{nn}}{\partial\omega}
         \right|_{\omega=\varepsilon_n^\mathrm{QP}}\right)^{-1} .

Between two neighbouring poles of :eq:`equ:gwbse:polesum` (for
:math:`\eta\to0`), :math:`\Sigma^\mathrm{c}` falls from :math:`+\infty` to
:math:`-\infty`, so :math:`f_n` decreases monotonically from
:math:`+\infty` to :math:`-\infty`: **there is exactly one root between any
two neighbouring poles**. A level therefore has as many roots as there are
poles, almost all of them with tiny weight -- satellites that sit just next
to a pole with a small residue. The root one wants is the *quasiparticle*
root with large :math:`Z` (typically 0.7–0.95 for valence levels). Away from
the gap, where the poles become dense and some carry large residues (above
roughly the lowest plasmon-like excitation energy of the system), a level can
have two or more roots of comparable weight. Then the QP picture of that
level is ambiguous, and which root is "the" QP energy is a matter of
definition (see :ref:`gwbse_satellites`).

The broadening :math:`\eta` makes :math:`f_n` smooth: each pole becomes a
steep but finite rise through zero with *positive* slope. These crossings are
not QP solutions (they have :math:`Z<0`) and are rejected.

``qp_solver``
=============

``grid`` (default)
    a search for sign changes of :math:`f_n` over a window, refinement of the
    roots found, and selection of one of them; described in the following;
``fixedpoint``
    Newton iteration from the current energy (``qp_sc_max_iter``,
    ``qp_sc_limit``), falling back to the grid search if it does not
    converge. It finds *a* root near the start, not necessarily the QP root,
    and is mainly useful for well-behaved levels.

If no root is found at all, the linearised solution
:math:`\varepsilon_n + Z_n f_n(\varepsilon_n)` is used and the level is
reported in the log.

The grid search
===============

For each level the search covers the window
:math:`[\omega_0-h, \omega_0+h]` around the current energy
:math:`\omega_0` (the DFT energy in the first iteration, the previous QP
input afterwards), with ``qp_full_window_half_width`` :math:`h` (default
0.75 Ha). It proceeds in three stages.

**Restricted window first** (``qp_restrict_search``, default ``true``). The
QP energy of an occupied level must lie below the Fermi level, that of a
virtual level above it. The first search therefore covers only the part of
the window below :math:`-`\ ``qp_zero_margin`` (occupied levels) or above
``qp_virtual_min_energy`` (virtual levels, default :math:`-0.1` Ha). If the
level itself lies outside these bounds -- for example an LUMO pulled below
:math:`-0.1` Ha by an environment -- the bound moves to the HOMO-LUMO
midpoint of the current energies. If the restricted search finds no
acceptable root, the full window is searched with the dense scan.

**Coarse adaptive scan** (``qp_grid_search_mode=adaptive_with_dense_fallback``,
default). Starting at the centre (in the first iteration the linearised
estimate, if it lies in the window), :math:`f_n` is evaluated at shells of
width ``qp_adaptive_shell_width`` (default 0.025 Ha) outwards on both sides,
all points in one batch (:ref:`gwbse_sigma`). In every interval with a sign
change, a local scan with 12 sub-steps locates the sign change(s) more
precisely; if there are several, the one nearest the centre is refined.

**Dense fallback.** If the adaptive scan yields no acceptable root, the
window is scanned uniformly with spacing ``qp_dense_spacing`` (default
0.002 Ha) and every sign change is refined. ``qp_grid_search_mode=dense``
uses only this scan, ``adaptive`` only the adaptive one.

Refinement and acceptance
=========================

Each sign-change interval is refined to a root by ``qp_root_finder``:

``brent`` (default)
    inverse quadratic interpolation with bisection as safeguard [Brent:1973]_;
    superlinear, a few evaluations per root;
``bisection``
    halving the interval; about :math:`\log_2(\Delta/\delta)` evaluations.

Both stop when :math:`|f_n| <` ``qp_sc_limit`` (default :math:`10^{-6}` Ha),
or when the interval is narrower than
:math:`\tfrac12 Z_\mathrm{min}\,` ``qp_sc_limit``. A root is **accepted** if

* :math:`|f_n(\omega)|\le` ``qp_sc_limit``,
* :math:`0.05 \le Z \le 1.5`.

The second stopping criterion guarantees that any root with an acceptable
weight also has an acceptable residual: near a root
:math:`|f_n'| = 1/Z \le 1/Z_\mathrm{min}`. (Earlier versions stopped at an
interval of ``qp_sc_limit`` itself; roots with :math:`Z` below about 0.5
could then end with a residual above the limit and were discarded although
they had been found -- with Brent, which returns an end of the interval
rather than its midpoint, more often than with bisection. In a CH4 test
case one member of a triply degenerate level ended on another root that
way.)

Choosing among several roots
============================

If several roots are accepted, one is chosen:

* by **weight** (first iteration, and always with
  ``qp_root_continuity=false``): the largest :math:`Z - 0.1\,|\omega-\omega_0|`,
  i.e. the largest weight with a slight preference for roots close to the
  current energy;
* by **continuity** (evGW iterations after the first, with
  ``qp_root_continuity=true``, the default): the root nearest to the
  current input energy :math:`\omega_0` of the level.

Continuity keeps every level on the branch it started on. This matters for
levels with two roots of comparable weight: choosing by weight anew in every
iteration lets such a level switch branches whenever a small change of
:math:`G` or :math:`W` tips the balance, and the self-consistency then never
settles (:ref:`gwbse_satellites`). For levels with one dominant root both
rules give the same result; on a QM9 molecule the results with and without
continuity agree in every printed digit.

Levels whose second-best accepted root has at least half the weight of the
chosen one are listed at the end of the GW step as *competing roots*
(chosen and alternative energy with their weights). These are the levels
whose QP energy depends on the choice of rule.

Accuracy of the roots and the evGW limit
========================================

The roots are located to within about :math:`Z\,` ``qp_sc_limit``. This
error enters the evGW residual (next section), so ``qp_sc_limit`` must be
clearly smaller than the evGW limit ``sc_limit``. With equal values (both
:math:`10^{-5}` Ha, the former defaults) the residual of the last iterations
is dominated by the rounding of the roots: it jumps between a few discrete
values around the limit, and whether a run "converges" depends on that
rounding rather than on the fixed point. The default ``qp_sc_limit`` is
therefore :math:`10^{-6}` Ha, a decade below ``sc_limit``. With Brent the
cost of the extra digit is small (QM9 molecule, exact :math:`\Sigma`: QP
step 6.5 s of 134 s).

Cost of the QP step
===================

Per level and iteration the grid search needs about 60 evaluations for the
coarse scan (in one batch) and a few per refined root. The scan dominates.
With the batched evaluation, the QP step is a minor part of the cost of an
evGW iteration with exact :math:`\Sigma`, and about 40 s per iteration for
the 995 levels of the D5M dimer with PPM. The debug log reports the number
of evaluations per iteration (``QP diagnostics``), the timing summary the
total time.


.. _gwbse_evgw:

Eigenvalue self-consistency (evGW)
**********************************

``mode=evGW`` (default) replaces the orbital energies in :math:`G` and
:math:`W` by the QP energies and iterates until they are consistent;
``mode=G0W0`` stops after the first QP solve. Only the energies are updated,
the orbitals stay the DFT ones (for orbital self-consistency see
:ref:`gwbse_qsgw`).

The fixed-point problem
=======================

Write :math:`\mathbf{E}` for the vector of energies of the GW window and
:math:`\mathbf{Q}(\mathbf{E}_G,\mathbf{E}_W)` for the QP energies obtained
with the energies :math:`\mathbf{E}_G` in :math:`G` and :math:`\mathbf{E}_W` in
:math:`W`. evGW looks for

.. math::

   \mathbf{E} = \mathbf{Q}(\mathbf{E},\mathbf{E}) .

Levels outside the GW window follow through ``out_of_window_shift``. An
iteration costs one QP solve per level and, if :math:`W` is rebuilt, one
screening build -- which is by far the more expensive part (PPM: two
dielectric matrices and the rotation of the integrals; exact: the RPA pair
eigenproblem).

Convergence criterion
=====================

An iteration is judged by its **QP residual**

.. math::

   r = \max_n \left| Q_n - E_n^\mathrm{in}\right| ,

the largest difference between the QP solutions and the energies they were
solved with, i.e. before any mixing. Earlier versions tested the change of
the *mixed* energies between iterations. With Anderson mixing that change can
become small while the energies are not a fixed point -- when the mixing
stagnates or a level alternates between two roots -- and runs then stopped
at a residual of up to :math:`10^{-3}` Ha. The residual test cannot be fooled
that way; on the other hand some runs now take a few more iterations.

At the info level (``--verbose1``) every iteration logs the residual and the
five levels with the largest residual as ``level: input -> solution
(residual)``. This line shows directly what holds a run up: levels that
creep towards their solution, levels that jump between two values, or core
levels that change after every rebuild of :math:`W`.

Mixing
======

The input of the next iteration is obtained from the history of inputs and
QP solutions by ``mixing_order``:

``0``
    plain update, :math:`\mathbf{E}^\mathrm{in}\leftarrow\mathbf{Q}`;
``1``
    linear mixing with ``mixing_alpha`` (default 0.7);
:math:`>1` (default 3)
    Anderson mixing [Anderson:1965]_ with a history of that many previous
    iterations.

Anderson mixing assumes that the inputs and outputs in its history come from
one smooth map. A short history recovers quickly when that is not the case
-- after a rebuild of :math:`W` (the history is restarted then, see below),
or when a level switches root. With a long history (the former default 20)
a few levels switching between nearby roots locked the mixing on a D5M dimer:
it returned the same input for several iterations although the QP solutions
were clearly elsewhere. On a QM9 molecule histories of 3, 5 and 20 give
identical results and iteration counts.

Adaptive screening update
=========================

``screening_update`` decides when :math:`W` is rebuilt:

``every``
    in every iteration (the classic scheme);
``adaptive`` (default)
    :math:`G` is iterated with :math:`W` fixed, and :math:`W` is rebuilt once
    the energies have settled for that :math:`W`.

With ``adaptive``, after each QP solve with residual :math:`r` the
**drift** :math:`d = \max_n |E_n - E_n^W|` of the energies from those
:math:`W` was built with is compared:

* :math:`r <` ``sc_limit`` and :math:`d <` ``sc_limit``: converged;
* :math:`r <` ``sc_limit`` but :math:`d \ge` ``sc_limit``: rebuild :math:`W`;
* :math:`r <` ``screening_update_ratio`` :math:`\times\, d` (default 0.1),
  or ``screening_update_max_inner`` (default 10) iterations with the same
  :math:`W`: rebuild :math:`W`;
* otherwise: next iteration with the same :math:`W`.

The convergence test requires both a small residual and a :math:`W` built
from the final energies, so the result is the same evGW fixed point as with
``every``. Only the path differs: :math:`G` converges against a frozen
:math:`W` first, and :math:`W` is rebuilt when :math:`G` has caught up. The
Anderson history is restarted at every rebuild, because the map changes with
:math:`W`.

Typical savings: QM9 molecule (exact :math:`\Sigma`) 10 → 6 screening
builds, GW time 192 s → 124 s; D5M dimer (PPM, 995 levels) 11 builds for 36
QP iterations, GW 3.4 h instead of 6 h. The summary line ``evGW: N QP
iterations, M screening builds (adaptive)`` reports both numbers.

Each rebuild of :math:`W` moves the energies by a fraction of their drift
since the previous build (about 10 % on QM9 molecules, about 30 % on the
dimer), so the number of builds grows with the logarithm of the initial
drift over ``sc_limit``. Levels that hardly influence :math:`W` but react to
it -- deep core levels in the GW window -- can dominate the last builds: on
the D5M dimer the last six of eleven builds were spent on five core levels at
:math:`-8.8` Ha that changed by :math:`10^{-4}` Ha after every rebuild, while
nothing near the gap moved any more. ``ignore_corelevels`` or ``qp_emin``
remove such levels from the GW window.

Root tracking (``qp_root_tracking``, optional)
==============================================

With ``qp_root_tracking=true`` the QP search of a level after the first
iteration covers only a window around its root of the previous iteration
(four times its last change, between 0.01 and 0.1 Ha) and keeps the nearest
accepted root, unless its weight has dropped below half of the previous one;
otherwise the full search is used. Once converged, one QP solve with the full
search checks the roots; if it chooses other roots for some levels, they are
logged and the run continues with the full search.

Tracking makes the QP step much cheaper (most levels need only a short local
scan) and suppresses root switching after rebuilds of :math:`W`. It has two
drawbacks. The root a level gets now depends on its previous root, not only
on the current energies, so the QP map changes from iteration to iteration
for levels that sometimes fall back to the full search; Anderson mixing can
then stall (seen on the D5M dimer). And where a level has two
self-consistent solutions, tracking may end on another one than the full
search would. Root continuity (the default) has neither problem: it uses the
full window in every iteration but chooses the root nearest the current
input, a function of the current energies only.


.. _gwbse_qsgw:

Quasiparticle self-consistent GW (QSGW)
***************************************

evGW updates only the energies; the orbitals stay those of the DFT start, and
so does the dependence of the result on the functional. QSGW
[vanSchilfgaarde:2006]_ also updates the orbitals: it replaces the
exchange-correlation potential by a static, Hermitian approximation of the
self-energy and iterates until the orbitals and energies of that
Hamiltonian reproduce themselves. ``do_qsgw=true`` runs QSGW after the
perturbative step (``mode=G0W0`` or ``evGW``), whose QP energies are the
starting point. The QSGW energies replace the perturbative ones, and the
BSE is set up in the QSGW orbitals.

QSGW is available for closed shells (RKS) with ``sigma_integrator=ppm`` or
``exact``. With ``cda`` it stops with an error: the quadrature becomes
ill-conditioned as the orbitals rotate.

The QSGW Hamiltonian
====================

In the DFT orbitals of the GW window the QSGW Hamiltonian is

.. math::
   :label: equ:gwbse:hqsgw

   H^\mathrm{QSGW} = \operatorname{diag}(\varepsilon^\mathrm{DFT})
     - V^\mathrm{xc} - a\,K[\rho^\mathrm{DFT}]
     + \Delta J
     + U\,\tilde\Sigma\,U^\mathsf{T} ,

with :math:`U` the current QP orbitals (columns, expanded in the DFT
orbitals) and :math:`\tilde\Sigma` the static self-energy in those orbitals.
Its eigenvalues are the QP energies, its eigenvectors the new QP orbitals.
The terms:

* :math:`\tilde\Sigma = \Sigma^\mathrm{x} + \tilde\Sigma^\mathrm{c}`
  (+ :math:`\Sigma^\mathrm{reac}` in a polarizable environment, see
  :doc:`embedded_gwbse`). :math:`\Sigma^\mathrm{x}` is the full exchange of
  the QP orbitals, with the frozen-core exchange if ``core_exchange=true``.
  :math:`\tilde\Sigma^\mathrm{c}` is the static correlation part, see
  :ref:`gwbse_qsgw_static`.
* **Hybrid starts.** :math:`V^\mathrm{xc}` holds only the semilocal part of
  a hybrid functional with exact-exchange fraction :math:`a`; the term
  :math:`-a K[\rho^\mathrm{DFT}]` removes the hybrid's exchange (including
  the frozen core) once, so that QSGW from PBE and from PBE0 solve the same
  equations. With a frozen core and ``core_exchange=false`` the core
  exchange would stay at :math:`a` times its value and the result would
  depend on the functional; the log warns about it. Range-separated hybrids
  are refused.
* **Hartree update** (``qsgw_hartree_update``, default true).
  :math:`H^\mathrm{KS}` holds the Hartree potential of the DFT density;
  :math:`\Delta J = J[\rho^\mathrm{QP}] - J[\rho^\mathrm{DFT}]` replaces it
  by that of the QP density of the active occupied levels (the frozen core
  does not change). Without it the occupied orbitals rotate in a fixed
  Hartree potential, which is not self-consistent.
* **Screened Hartree update** (``qsgw_hartree_screened``, default true). In
  a polarizable environment :math:`H^\mathrm{KS}` also holds the potential
  of the dipoles the DFT density induced. The update then uses
  :math:`u = v + v^\mathrm{reac}`, the interaction :math:`W` is built
  from, so the induced dipoles follow the QP density. ``false`` keeps the
  environment's response to the DFT density, a control: with it the fixed
  point still depends on the start (6--40 meV median spread for methane in
  five Thole sites, depending on how far the starts are apart).

The loop
========

Each iteration

#. rotates the three-centre integrals to the current QP orbitals (both
   orbital indices) and rebuilds :math:`W` with the current QP energies;
#. evaluates :math:`\tilde\Sigma` and :math:`\Delta J` in these orbitals
   and transforms them to the DFT orbitals;
#. mixes :math:`\tilde\Sigma` (the whole matrix) with the history of
   previous iterations: Anderson mixing over ``mixing_order`` previous
   iterations followed by linear mixing with ``mixing_alpha``;
   ``mixing_order=0`` is linear mixing alone;
#. diagonalises :eq:`equ:gwbse:hqsgw` for the next QP energies and
   orbitals.

The loop has converged when the self-energy it produces differs from the one
it was given by less than ``qsgw_sc_limit`` (largest matrix element) and the
QP energies change by less than ``qsgw_sc_limit``; the result is then the
eigensystem of the produced Hamiltonian. After ``qsgw_max_iterations``
without convergence the log warns and the last iterate is kept.

Every iteration logs ``max|dE_QP|`` and ``max|Sigma_out - Sigma_in|``. With
the Hartree update or an environment a second line gives the largest change
of the Hartree potential (bare and environment parts), the change of
:math:`\Sigma^\mathrm{reac}` for HOMO and LUMO, and the largest rotation
between the occupied DFT and QP spaces in degrees.

**Virtual levels far from their DFT energy.** High virtual levels often have
perturbative corrections of several eV and destabilise the loop. Before the
loop, the GW window is cut below the first virtual level whose perturbative
correction exceeds ``qsgw_max_virt_correction`` (default 0.5 Ha); the
levels above keep their perturbative energies and DFT orbitals. The log
reports the cut.

.. _gwbse_qsgw_static:

The static self-energy and its regularization
=============================================

:math:`\tilde\Sigma^\mathrm{c}` has to be static and Hermitian, but the
correlation self-energy is neither. Off the diagonal it generalises the pole
sum :eq:`equ:gwbse:polesum` to

.. math::

   \Sigma^\mathrm{c}_{nn'}(\omega) = \sum_{p}
     \frac{b_{nn'p}\,(\omega - z_{p})}{(\omega - z_{p})^2 + \eta^2},
   \qquad b_{nn'p} = c_p\, r_{np}\, r_{n'p} ,

with :math:`p = (m,s)` running over the poles, :math:`r_{np}` the residues
and :math:`c_p` the mode weights of the self-energy model. Two static forms
are implemented (``qsgw_regularization``).

``eta`` (default)
    the symmetrised self-energy at the QP energies,

    .. math::

       \tilde\Sigma^\mathrm{c}_{nn'} = \tfrac12\left[
         \Sigma^\mathrm{c}_{nn'}(E_n) + \Sigma^\mathrm{c}_{nn'}(E_{n'})\right] ,

    with the broadening ``eta``. A pole :math:`z_p` close to a QP energy
    gives a term up to :math:`1/(2\eta)` in size, so with the default
    :math:`\eta = 10^{-3}` Ha the Hamiltonian changes abruptly as levels
    move past poles.

``srg``
    the form of the similarity-renormalization-group derivation of Marie
    and Loos [Marie:2023]_ (their eq. 44): with
    :math:`\Delta_{np} = E_n - z_p`,

    .. math::

       \tilde\Sigma^\mathrm{c}_{nn'} = \sum_p b_{nn'p}\,
         \frac{\Delta_{np} + \Delta_{n'p}}{\Delta_{np}^2 + \Delta_{n'p}^2}
         \left[1 - e^{-s(\Delta_{np}^2 + \Delta_{n'p}^2)}\right] ,

    on the diagonal
    :math:`b_{nnp}\,(1 - e^{-2s\Delta_{np}^2})/\Delta_{np}`. It is
    Hermitian without symmetrisation, and every term stays finite: a pole at
    distance :math:`\Delta` from both levels contributes at most
    :math:`\sim\sqrt{s}`, and terms with :math:`|\Delta| \gg s^{-1/2}`
    are the undamped :math:`1/\Delta`. The flow parameter
    ``qsgw_srg_s`` (in Ha\ :sup:`-2`, default 500) plays the role of
    :math:`1/(2\eta^2)`: :math:`s = 200` corresponds to
    :math:`\eta = 0.05` Ha, the default 500 to :math:`\eta \approx 0.032` Ha.
    For :math:`s \to \infty` the diagonal becomes that of the ``eta`` form at
    :math:`\eta \to 0`; for :math:`s \to 0` the correlation part vanishes
    (Hartree-Fock). Marie and Loos recommend :math:`s` = 500 or 1000; in their
    molecular test set qsGW did not converge for larger :math:`s`.

**Why the regularization matters.** QSGW is a nonlinear fixed-point problem
and can have several solutions. With small :math:`\eta` it often does: for
methane/3-21G with exact :math:`\Sigma` and all levels in the window, starts
from 20 randomly rotated orbital sets reach another fixed point in 3--10
cases, with HOMO or LUMO up to 0.18 eV apart; 10--25 % of the runs do not
converge in 100 iterations, and the others need 45--90. A start-dependent result is then not
a property of the method but of the solution the run happened to find. Both
regularizations remove this, at a price that differs:

.. list-table:: QSGW of methane/3-21G in vacuum, exact :math:`\Sigma` (eV)
   :header-rows: 1
   :widths: 30 20 20 30

   * - Setting
     - HOMO
     - Gap
     - Unique fixed point
   * - ``eta``, :math:`\eta = 10^{-3}` Ha
     - -17.6708
     - 24.9863
     - no
   * - ``eta``, :math:`\eta = 0.05` Ha
     - -17.6623
     - 24.9711
     - yes
   * - ``srg``, :math:`s = 100`
     - -17.6654
     - 24.9814
     - yes
   * - ``srg``, :math:`s = 20000`
     - -17.6654
     - 24.9819
     - yes

With ``eta``, a fixed point is unique only once :math:`\eta` is large, and
then the QP energies depend on :math:`\eta` (here 15 meV in the gap between
:math:`10^{-3}` and 0.05 Ha). With ``srg``, the energies change by less than
1 meV between :math:`s` = 100 and 20000, every converged start reaches the
same fixed point for :math:`s` = 200--1000 (in vacuum and in an
environment), and the loop needs 10--50 iterations, fewer for smaller
:math:`s`. For water (def2-TZVP, exact :math:`\Sigma`) the picture is the
same: PBE and PBE0 starts agree within 0.7 meV for :math:`s` = 200--1000,
against 32 meV at :math:`\eta = 0.01` Ha.

**SRG and the** ``eta`` **form are different methods.** Off the diagonal the
two static operators differ even where no regularization acts: for a pole
far from both levels, ``srg`` gives
:math:`(\Delta_{np} + \Delta_{n'p})/(\Delta_{np}^2 + \Delta_{n'p}^2)`, the
``eta`` form :math:`\tfrac12(1/\Delta_{np} + 1/\Delta_{n'p})`. The QP energies
therefore differ, by an amount that depends on the system: 5 meV in the HOMO
of methane, but 0.18 eV for water (HOMO :math:`-12.71` eV with ``srg`` for
any :math:`s` from 50 to 20000, :math:`-12.89` to :math:`-12.91` eV with
``eta`` for :math:`\eta` = 0.005--0.05 Ha). Marie and Loos find the same
direction over 50 molecules: SRG-qsGW ionization potentials are 0.06 eV
lower on average than qsGW ones (water: 0.10 eV) and closer to
:math:`\Delta`\ CCSD(T). Do not mix results of the two forms in one
comparison.

The disabled unit tests ``qsgw_fixed_point_seed_scan`` and
``qsgw_regularization_scan`` (``test_gw``) reproduce these numbers; run
them with ``--run_test=gw_test/<name>``.

Recommended settings
====================

* Use ``qsgw_regularization=srg`` with the default ``qsgw_srg_s=500``, and
  report :math:`s` and the form. To quote a regularization uncertainty,
  repeat with :math:`s` = 200 and 1000. Much larger :math:`s` brings back
  near-singular terms: for water the loop needs 74--172 iterations at
  :math:`s = 20000`, and the two starts differ by 1.3--1.5 meV for
  :math:`s` = 5000--20000; Marie and Loos
  report failed convergence for 11 of 50 molecules at :math:`s = 5000`.
* Check that the fixed point does not depend on the start before
  interpreting a starting-point dependence -- for instance QSGW from PBE and
  from PBE0 should agree to within ``qsgw_sc_limit``.
* Keep ``core_exchange=true`` with a frozen core, in particular from hybrid
  starts.
* With environment screening keep ``qsgw_hartree_screened=true``.
* Keep the default mixing (``mixing_order=3``, ``mixing_alpha=0.7``). A long
  Anderson history stalls: for water with ``srg``, order 20 with
  :math:`\alpha = 0.2` left a residual of 0.3 eV after 200 iterations where
  the default converges in 20.
* If the loop does not converge, check the virtual cut
  (``qsgw_max_virt_correction``) and the mixing first; ``mixing_order=0``
  with a small ``mixing_alpha`` shows whether the problem is the mixing or
  the fixed point.

.. list-table:: QSGW options (``gwbse.gw``)
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Default
     - Meaning
   * - ``do_qsgw``
     - false
     - run QSGW after G0W0 or evGW
   * - ``qsgw_regularization``
     - ``eta``
     - static self-energy: ``eta`` or ``srg``
   * - ``qsgw_srg_s``
     - 500
     - SRG flow parameter (Ha\ :sup:`-2`)
   * - ``eta``
     - 1e-3
     - broadening of the ``eta`` form (Ha); also used by the perturbative
       step
   * - ``qsgw_hartree_update``
     - true
     - Hartree potential of the QP density
   * - ``qsgw_hartree_screened``
     - true
     - with environment screening: update with :math:`v + v^\mathrm{reac}`
   * - ``qsgw_sc_limit``
     - 1e-5
     - convergence threshold (Ha)
   * - ``qsgw_max_iterations``
     - 50
     - maximum number of iterations
   * - ``qsgw_max_virt_correction``
     - 0.5
     - cut of the window above the first virtual level with a larger
       perturbative correction (Ha)
   * - ``mixing_order``, ``mixing_alpha``
     - 3, 0.7
     - mixing of :math:`\tilde\Sigma` (shared with evGW)


.. _gwbse_satellites:

Levels in the satellite region
******************************

Far enough above the gap (or below it, for deep occupied levels), the poles
of :eq:`equ:gwbse:polesum` become dense and some carry large residues. A
level there can have two or more roots of comparable weight, and the QP
picture of that level breaks down; multiple solutions of the QP equation and
the discontinuities they cause are well documented [vanSetten:2015]_
[Veril:2018]_. Three things can then happen in evGW.

**Root switching.** If the root is chosen by weight in every iteration, a
level with two roots of similar :math:`Z` takes whichever is ahead at the
moment; small changes of :math:`G` or, after a rebuild, of :math:`W` move it
back and forth. The residual then stays at the distance between the roots.
Root continuity (default) prevents this.

**Several fixed points.** The two roots can each belong to a self-consistent
solution: the level is consistent on either branch. Which one a run ends on
depends on the path -- the root rule, the mixing, even rounding. An example
is an open-shell cation in a small basis whose beta level 12 has
self-consistent solutions at :math:`+2.41` and :math:`+2.73` Ha; different
settings reached either. Both are valid evGW solutions; the QP energy of
such a level is not well defined.

**No fixed point.** For some levels the root near the input moves away as
the input approaches it, because the level's own energy shifts the poles it
sees. A performance run on C60 with the generous default window showed
degenerate sets of virtual levels about 24 eV above the LUMO alternating
between roots 25–29 mHa apart: from input 0.780 Ha the solution was 0.808 Ha,
from input 0.808 Ha it was 0.780 Ha. No setting of the root search can make
such a level converge, and the run stops after ``sc_max_iter`` iterations
with a warning.

The levels concerned are far from the gap. They hardly influence the
low-lying excitations, but evGW requires every level of the window to
converge, and they cost time in every iteration. The remedy is to keep them
out of the GW window (``qp_emax``, see :ref:`gwbse_windows`); above the
window their energies follow ``out_of_window_shift``. The top-5 residual
lines and the competing-roots list identify them.


.. _gwbse_bse:

Bethe-Salpeter equation
***********************

Hamiltonian and solvers
=======================

The BSE is built in the product space of the BSE window, with the QP
energies of the GW window (diagonal unless ``use_Hqp_offdiag``), the exchange
term with the bare Coulomb interaction and the direct term with the static
screened interaction :math:`W(\omega=0)`. ``exctotal`` sets the number of
states.

With ``useTDA=true`` only the resonant block is solved (Tamm-Dancoff
approximation). Otherwise the full BSE, whose eigenvalue problem for real
orbitals is equivalent to the symmetric problem

.. math::

   (A-B)^{1/2}(A+B)(A-B)^{1/2}\,\mathbf{T} = \Omega^2\,\mathbf{T} ,

is solved in this Hermitian form by a dedicated Davidson solver; it never
forms the dense blocks for the exchange and QP terms. ``davidson`` sets the
correction method, tolerance, subspace update and iteration limit. The log
prints every Davidson iteration with the largest residual and the number of
converged roots; the solver stops once the requested roots are converged.

Caches for the direct term
==========================

The direct term needs :math:`W(0)` between electron-electron and hole-hole
products. It can be applied from the three-centre integrals in every Davidson
iteration or from a precomputed matrix. The precomputed matrices are kept as
lower-triangle tiles (half the memory and half the flops to build) and
shared between the solves and the analysis of the states. Whether they
are kept depends on the memory budget (``--memory`` in GB on the command
line); the log reports the decision, e.g. ``BSE direct terms (lower
triangle, 0.0171 GB each): Hd kept; Hd2 kept``.

Dynamical screening
===================

``dyn_screen_max_iter`` > 0 adds a perturbative correction for the frequency
dependence of :math:`W` to each BSE state, iterated to ``dyn_screen_tol``.
The corrected energies are printed after the static ones (in the order of the
static roots).


Performance and memory
**********************

Where the time goes
===================

The log ends the GW step with a timing summary. Typical shares:

* **PPM, large system** (D5M dimer, 995 GW levels, 11 builds): rotation of
  the integrals into the PPM frame 59 %, dielectric matrices 26 %, QP
  equations 13 %; off-diagonal :math:`\Sigma` another 13 minutes after GW.
* **Exact** :math:`\Sigma` **, small molecule** (QM9, def2-TZVP): the RPA
  pair eigenproblem 88 %, residues 5 %, QP equations 6 %.

The number of screening builds is the main lever; the adaptive update,
the default settings and a GW window without satellite-region and deep-core
levels reduce it.

For many small molecules (data sets), the dense eigensolvers run far below
peak on many threads; running several molecules in parallel with fewer
threads each gives more molecules per hour than one molecule on all threads.

Memory
======

``--memory`` (GB) is the budget for choices that trade memory for speed
(BSE caches, batch sizes); it does not limit what a calculation needs
anyway. The log reports the resident and peak memory after the three-centre
integrals and around the first screening build. The three-centre integrals
of the GW window dominate for large systems (D5M dimer: about 200 GB).
Exact :math:`\Sigma` needs the RPA pair matrix
(:math:`(N_\mathrm{occ}N_\mathrm{virt})^2`) and the residues.

Platforms
=========

The dense linear algebra uses MKL where available (``dsyevd``, threaded
``dsyrk``), Apple Accelerate's LAPACK on macOS 13.3 and later, and Eigen
otherwise; the log names the backend. CUDA support has been removed.


Troubleshooting
***************

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Symptom in the log
     - Likely cause
     - What to do
   * - ``WARNING! GW-self-consistency cycle not converged``; the top-5
       residual lines show the same few levels jumping by tens of mHa
     - levels in the satellite region with two roots (switching, or no
       fixed point)
     - cap the GW window (``qp_emax``) below these levels; keep
       ``qp_root_continuity=true``
   * - residual stuck at one value for many iterations, inputs of the listed
       levels repeat exactly
     - Anderson mixing stalled by a discontinuous QP map
     - shorter ``mixing_order`` (default 3); avoid ``qp_root_tracking``
   * - last screening builds only move core levels (top-5 lines list levels
       far below the HOMO)
     - core levels in the GW window react to every new :math:`W`
     - ``ignore_corelevels=RPA`` (with ``core_exchange``) or ``qp_emin``
   * - residual of the last iterations alternates between a few discrete
       values near ``sc_limit``
     - roots located no more accurately than the evGW limit
     - ``qp_sc_limit`` at least ten times smaller than ``sc_limit``
       (default)
   * - ``Competing roots`` listed for a level of interest
     - two roots of comparable weight; the QP energy depends on the
       choice rule
     - inspect :math:`\Sigma(\omega)` (``sigma_plot``); compare runs with
       ``qp_root_continuity`` true and false
   * - ``Restricted QP search failed ... forcing full dense window``
     - the level has no acceptable root on its side of the Fermi level
       within the window
     - usually harmless; frequent occurrences point to the satellite region
       or to a level pulled across the gap (environment)
   * - frozen-core and all-electron results differ by tenths of an eV
     - core exchange missing (old versions, or ``core_exchange=false``)
     - ``core_exchange=true`` (default)
   * - PPM log reports modes with negative squared frequency
     - the PPM fit at the imaginary frequency failed for these modes
     - harmless if few and with small weight; compare with exact
       :math:`\Sigma` for a small system if in doubt

Reproducing results of older versions
=====================================

The defaults changed in this version. To reproduce earlier numbers, set the
previous values explicitly::

    <gw>
      <screening_update>every</screening_update>
      <mixing_order>20</mixing_order>
      <qp_root_continuity>false</qp_root_continuity>
      <qp_root_finder>bisection</qp_root_finder>
      <qp_sc_limit>1e-5</qp_sc_limit>
      <rebuild_3c_freq>5</rebuild_3c_freq>
      <core_exchange>false</core_exchange>
    </gw>

Even then results can differ where the old convergence test stopped early,
where the old root refinement discarded small-weight roots, and for frozen
core (the old frozen-core Σx was incomplete).


Reference: evGW and QP settings
*******************************

All in the ``gw`` subsection unless noted.

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Option
     - Default
     - Meaning
   * - ``mode``
     - ``evGW``
     - ``G0W0`` or ``evGW``
   * - ``sigma_integrator``
     - ``ppm``
     - ``ppm``, ``exact`` or ``cda``
   * - ``eta``
     - 1e-3
     - broadening of the poles (Ha)
   * - ``sc_limit``
     - 1e-5
     - evGW limit on the QP residual (and on the drift with ``adaptive``), Ha
   * - ``sc_max_iter``
     - 50
     - maximum number of QP iterations
   * - ``screening_update``
     - ``adaptive``
     - when :math:`W` is rebuilt: ``every`` or ``adaptive``
   * - ``screening_update_ratio``
     - 0.1
     - rebuild once the residual is below this fraction of the drift
   * - ``screening_update_max_inner``
     - 10
     - rebuild at the latest after this many iterations with the same
       :math:`W`
   * - ``mixing_order``
     - 3
     - 0 plain, 1 linear, >1 Anderson history
   * - ``mixing_alpha``
     - 0.7
     - mixing parameter
   * - ``qp_solver``
     - ``grid``
     - ``grid`` or ``fixedpoint``
   * - ``qp_grid_search_mode``
     - ``adaptive_with_dense_fallback``
     - ``adaptive``, ``dense`` or both
   * - ``qp_full_window_half_width``
     - 0.75
     - half-width of the search window (Ha)
   * - ``qp_adaptive_shell_width``
     - 0.025
     - spacing of the coarse scan (Ha)
   * - ``qp_dense_spacing``
     - 0.002
     - spacing of the dense scan (Ha)
   * - ``qp_restrict_search``
     - true
     - search occupied levels below and virtual levels above the Fermi
       level first
   * - ``qp_root_finder``
     - ``brent``
     - ``brent`` or ``bisection``
   * - ``qp_sc_limit``
     - 1e-6
     - residual of an accepted root (Ha)
   * - ``qp_root_continuity``
     - true
     - after the first iteration choose the root nearest the input
   * - ``qp_root_tracking``
     - false
     - local search near the previous root, verified at the end
   * - ``out_of_window_shift``
     - ``max``
     - shift of levels outside the GW window: ``max``, ``nearest``,
       ``linear``
   * - ``core_exchange``
     - true
     - frozen core: exchange with the core levels in :math:`\Sigma_x`
   * - ``rebuild_3c_freq``
     - 0
     - rebuild the three-centre integrals every this many iterations
       (0: never)
   * - ``qp_emax``, ``qp_emin``, ``bse_emax``, ``bse_emin`` (``gwbse``)
     - 0
     - energy caps of the GW and BSE windows (eV, 0: off)
   * - ``ignore_corelevels`` (``gwbse``)
     - ``none``
     - ``RPA``, ``GW``, ``BSE`` or ``none``
