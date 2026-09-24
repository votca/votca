Polarizable Ewald Embedding
###########################

A charge carrier in a molecular solid does not sit in vacuum. Its energy is
shifted by the electrostatic potential of every other molecule in the
sample, and by the polarization those molecules undergo in response to the
carrier's own field. Both contributions are long ranged: the electrostatic
one falls off only as :math:`1/r` for a charged carrier, and the
polarization response is itself driven by that slowly decaying field. A
calculation that simply truncates the environment at some radius does not
converge to the right answer, it converges to an answer that depends on
where the truncation was placed.

The machinery described here computes those shifts for a periodic sample,
by combining an Ewald lattice sum over a self-consistently polarized
background with an explicit, quantum-mechanical treatment of the molecule
carrying the charge. It is used in two ways:

* as a purely classical calculation, in which every molecule including the
  central one is represented by distributed multipoles, and

* as a three-region QM/MM calculation, in which the central molecule is
  described by DFT, a shell of neighbours is treated as explicitly
  polarizable point multipoles, and everything beyond that is the periodic
  Ewald background.

The two share all of their machinery below the level of the central
molecule's description, which is deliberate: it means the classical
calculation is a reference for the QM/MM one, and a disagreement between
them is a statement about the description of that molecule rather than
about the environment.


The Ewald decomposition
***********************

The quantity to be summed is the electrostatic interaction between a set of
multipoles in a periodic cell and all of their periodic images,

.. math::
    :label: equ:ewald:latticesum

    E = \frac{1}{2}\sum_{\mathbf{T}}\ {\sum_{i,j}}'
        Q_i \, \hat{T}(\mathbf{r}_i - \mathbf{r}_j - \mathbf{T}) \, Q_j ,

where :math:`\mathbf{T}` runs over lattice translations, :math:`Q_i`
denotes the multipole moments of site :math:`i`, :math:`\hat{T}` is the
interaction tensor, and the prime excludes :math:`i = j` in the
:math:`\mathbf{T} = 0` cell. This sum is only *conditionally* convergent:
its value depends on the order in which the terms are taken, which is
another way of saying it depends on the shape of the macroscopic sample and
on the boundary condition applied at its surface.

The Ewald construction splits the :math:`1/r` kernel with the error
function,

.. math::
    :label: equ:ewald:split

    \frac{1}{r} = \underbrace{\frac{\mathrm{erfc}(\alpha r)}{r}}_{\text{short ranged}}
                + \underbrace{\frac{\mathrm{erf}(\alpha r)}{r}}_{\text{smooth}} ,

and evaluates the two pieces in the spaces where each converges quickly.
The screened part is summed directly over neighbouring images in real
space; the smooth part is summed in reciprocal space, where its Fourier
transform decays as :math:`\exp(-k^2/4\alpha^2)`. The result is assembled
from four contributions:

Real space
    The erfc-screened interaction, summed over all pairs and lattice
    translations within a distance cutoff. Implemented by
    ``EwaldRealSpaceSum``.

Reciprocal space
    The erf-screened interaction, as a sum over reciprocal lattice vectors
    weighted by :math:`\exp(-k^2/4\alpha^2)/k^2`. Implemented by
    ``EwaldReciprocalSpaceSum``, which builds a structure factor once per
    call and replays it against every target. The :math:`k = 0` term is
    omitted; see `The k = 0 term and the potential's gauge`_.

Shape (surface) term
    The uniform depolarizing field that encodes the boundary condition of
    :eq:`equ:ewald:latticesum`, expressed through the total dipole moment
    :math:`\mathbf{M} = \sum_j (q_j \mathbf{r}_j + \boldsymbol{\mu}_j)` of
    the cell [DeLeeuw:1980]_. Under vacuum boundary conditions,

    .. math::

        \mathbf{E} = -\frac{4\pi}{3V}\mathbf{M} \quad\text{(cube or sphere)},
        \qquad
        \mathbf{E} = -\frac{4\pi}{V}M_z\,\hat{\mathbf{z}} \quad\text{(slab)} .

    Implemented by ``EwaldShapeCorrection``. This is not a pairwise sum: it
    is a mean-field property of the whole sample's surface polarization
    acting back on itself, so it includes the target's own segment.

Self-interaction removal
    The reciprocal sum, being a sum over the full periodic density, includes
    each site's interaction with itself. For a point charge this term
    vanishes by symmetry. For a *static dipole* it does not: a dipole's own
    erf-screened field at its own position is finite and equals
    :math:`\tfrac{4}{3}\alpha^3/\sqrt{\pi}`. Omitting it makes the total
    depend on :math:`\alpha`, and this was a real defect in an earlier
    version of this code.

.. _ewald-alpha-invariant:

The :math:`\alpha` invariant
============================

The splitting parameter :math:`\alpha` in :eq:`equ:ewald:split` is
arbitrary. It controls how work is divided between the real- and
reciprocal-space sums, and nothing else. **The total must therefore be
independent of it**, and any dependence is a bug rather than a tolerance.

This is the single most useful diagnostic in the whole implementation,
because it is sensitive to a *missing term* rather than to a wrong factor.
A term that is absent from one sum but present in the other will leave a
residue that grows or shrinks as the split moves, even when every
individual number looks plausible. Both defects mentioned above — the
dipole self-term, and an omitted erf correction — were found this way, and
the current implementation is flat to the tenth digit across a factor of
two in :math:`\alpha` while the individual channels move by factors of
twenty to sixty.

Anyone changing this code should run an :math:`\alpha` scan before
believing a result. Convergence parameters must be held *fixed* while
:math:`\alpha` varies, not scaled with it: a truncation error that tracks
:math:`\alpha` imitates exactly the failure the test is looking for.


Segments, fragments and charge states
*************************************

The classical representation of a molecule is a ``PolarSegment``: a list of
``PolarSite`` objects, each carrying a position, a set of permanent
multipole moments up to the requested rank, a polarizability tensor, and an
induced dipole. Segments are built from the mapping file by the polar
mapper, which associates each coarse-grained segment with a ``.mps`` file
per charge state.

Multipoles are read from ``.mps`` files, in the format used by GDMA and
related tools. Two conventions in that format are worth stating explicitly
because they are easy to get wrong:

* the dipole line is ordered **z x y**, not x y z;
* multipoles are in :math:`e\,a_0` and polarizabilities in
  :math:`\mathrm{\AA}^3`, while the code works internally in atomic units
  throughout (bohr, Hartree).

Each segment is registered for the charge states it is needed in —
``Neutral``, ``Electron``, ``Hole`` — and these are held in an
``EwaldRegistry``, keyed by segment id and state. The registry is the single
source of truth for the periodic density: it replaces the legacy design in
which several containers shared raw pointers to the same segments and their
charge state depended on which container had most recently been asked.

A site's permanent and induced moments are kept separate throughout, and
they enter the energy through different channels. This matters for the
bookkeeping described in `Energy accounting`_, and for the Thole damping
described next, which applies to induced-dipole interactions and not to
permanent ones.


The polarizable background
**************************

Before any job is run, the periodic environment must be polarized
self-consistently. This is the job of the ``ewaldbackground`` calculator,
which is run once per frame:

.. code-block:: bash

    xtp_run -e ewaldbackground -o ewaldbackground.xml -f state.hdf5

Each site :math:`i` carries an induced dipole
:math:`\boldsymbol{\mu}_i` responding to the total field at its position,

.. math::
    :label: equ:ewald:induction

    \boldsymbol{\mu}_i = \boldsymbol{\alpha}_i \Big(
        \mathbf{F}^{\text{perm}}_i
        + \sum_{j \neq i} \mathbf{T}^{\text{Thole}}_{ij}\,\boldsymbol{\mu}_j
      \Big) ,

which is a linear system in the induced dipoles, solved to
self-consistency. The field :math:`\mathbf{F}^{\text{perm}}` is the full
periodic Ewald field of the permanent multipoles; the coupling
:math:`\mathbf{T}^{\text{Thole}}` between induced dipoles is damped at
short range by the exponential Thole model [Thole:1981]_, which prevents
the polarization catastrophe that an undamped point-dipole model suffers
when two polarizable sites approach each other.

The damping factors depend on both sites' polarizabilities through
:math:`u^3 = a\,r^3 \sqrt{\alpha_i \alpha_j}^{-1}`, so a site with no
polarizability is read as *complete* overlap and damped maximally. That is
the correct limit for two real sites and the wrong one for a field
evaluation point; see `Coupling to the QM region`_.

Two solvers are available. The default is a preconditioned conjugate
gradient; a Jacobi over-relaxation scheme reproducing the legacy SOR
iteration is retained for direct comparison. The PCG solver additionally
reports a Lanczos estimate of the smallest eigenvalue of the interaction
operator, which is a direct algebraic certificate of whether the system is
positive definite rather than an inference from residual behaviour.

The result is written to a checkpoint, by default
``ewaldbackground.hdf5``, containing

* the ``EwaldRegistry`` — every segment, in every registered charge state,
  with its converged induced dipoles, and
* an ``ewald_parameters`` group holding
  :math:`\alpha`, :math:`k_\text{max}`, :math:`r_\text{min}`, the field
  tolerance, the Thole parameter :math:`a`, the screening factor, the
  shape, and the simulation box.

**The parameters are part of the checkpoint on purpose.** See
`Parameter inheritance`_.


The three-region QM/MM setup
****************************

A QM/MM job is built from three regions, which must be enumerated in this
order:

.. list-table::
   :header-rows: 1
   :widths: 8 22 70

   * - id
     - type
     - role
   * - 0
     - ``qmregion``
     - The molecule carrying the charge, described by DFT. Holds the
       segment whose site energy is wanted, in the charge state the job
       specifies.
   * - 1
     - ``polarregion``
     - An explicit shell of neighbours, as polarizable point multipoles.
       These respond to the QM density and to the background, and they
       polarize each other.
   * - 2
     - ``ewaldregion``
     - The periodic background, read from the ``ewaldbackground``
       checkpoint. **Frozen**: nothing in the job repolarizes it.

An ``ewaldregion`` may not be region 0. It owns no segments carved out of
the job's topology — it represents the whole periodic cell — so the
recentring that ``JobTopology`` applies has no segment of its own to key
on.

Because the background is frozen, ``EwaldRegion`` is always ``Converged()``
and its ``Reset()`` is a no-op. Including it can never prevent the
inter-region SCF loop from terminating. Its ``Interactwith*`` methods return
zero, and this is the correct physics rather than a stub: no other region
polarizes it. The influence in the other direction is implemented by *that*
region's own ``InteractwithEwaldRegion``, following the receiver-pull
convention used throughout ``Region::ApplyInfluenceOfOtherRegions``.

A two-region ``qmregion`` + ``ewaldregion`` job is also valid and is useful
as a diagnostic, since it removes the classical coupling entirely.

.. _ewald-foreground:

Foreground declaration and suppression
======================================

The segments treated explicitly — by the QM region and by the polar region
together — are also present in the background, since the background is the
whole cell. They must not be counted twice. Removing them is not simply a
matter of dropping them from the sums:

* in **real space** the coincident copy of each explicit segment is
  suppressed, while its periodic *images* are kept. Deleting the images too
  would replace one carved-out cavity with a lattice of vacancies.
* in **reciprocal space** nothing is held out. A :math:`k`-space sum runs
  over the whole periodic density and cannot have a hole cut in it.
  Instead, the *erf-screened* energy of exactly those coincident copies is
  subtracted afterwards, which removes precisely what the reciprocal sum
  put back.

The copies to be suppressed are the union over *all* regions that own
segments, and no single region knows that union. ``JobTopology`` therefore
calls ``EwaldRegion::RegisterForeground`` once, after building every region
and before evaluating any of them. Without this, a QM/MM job would leave the
QM segment's neutral background copy sitting underneath the QM density — a
ghost molecule in every sum, with no symptom to notice it by. The suppressed
count is reported per job and checked against the expected total.

The erf correction uses the background's *own* multipoles and induced
dipoles, not the job's charge state, because that is what the reciprocal sum
actually placed there. Using the job's state would remove something that was
never added.


Energy accounting
*****************

Each pairwise contribution is computed exactly once, and which region books
it follows a single rule: **the region further inside owns the energy**. A
site accumulates the field from regions further out in ``V()``, which
carries energy, and from regions further in in ``V_noE()``, which does not.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - contribution
     - where it is booked
   * - QM internal
     - the DFT total energy
   * - QM ↔ polar
     - the DFT total, via the external multipole matrix (electrons) and
       ``ExternalRepulsion`` (nuclei). ``PolarRegion::InteractwithQMRegion``
       returns zero by this convention.
   * - QM ↔ background
     - the DFT total, via the potential on the integration grid (electrons)
       and a scalar added to :math:`E_0` (nuclei)
   * - polar internal (permanent)
     - ``E_static_static``
   * - polar ↔ background, permanent × permanent and permanent × induced
     - ``E_static_ext``, returned by ``EwaldRegion::ApplyFieldTo``
   * - polar induced × background
     - ``E_polar_ext``, as :math:`\sum \boldsymbol{\mu}\cdot\mathbf{V}`
   * - polar induction self-energy
     - ``E_polar_internal``
   * - background internal
     - **not computed** — see below

The background's own internal energy is deliberately absent. It is a
constant, independent of the job's charge state, so it cancels exactly in
any charge-state difference. Site energies are therefore correct; an
absolute total energy is not, and should not be quoted as one.

One sign convention deserves attention. The Ewald code and the region
framework store :math:`\mathbf{V}` with *opposite* signs — the background
solver builds its right-hand side as :math:`+\mathbf{V}`, the polar region
as :math:`-(\mathbf{V} + \mathbf{V}_{\text{noE}})`. ``ApplyFieldTo``
negates its own contribution at the boundary so that the polar region's
induction runs in the right direction. Both conventions are internally
consistent; only the boundary between them needs care.

Site energies
=============

The quantity of interest is the difference between charge states, referenced
to the same molecule in vacuum:

.. math::
    :label: equ:ewald:siteenergy

    \Delta E_{\text{env}} =
      \big[E(h) - E(n)\big]_{\text{embedded}} -
      \big[E(h) - E(n)\big]_{\text{vacuum}} .

Using the same geometries in both brackets makes the intramolecular
relaxation cancel exactly, leaving the environment's response to the charge
at its own geometry. In practice this means running the same job file
against a ``qmregion``-only setup to obtain the vacuum reference.


Coupling to the QM region
*************************

The background reaches the QM region as a *potential*, evaluated on the
same numerical integration grid the DFT calculation already uses:

.. math::

    E_{\text{QM-bg}} =
      -\int \rho(\mathbf{r})\,\phi(\mathbf{r})\,\mathrm{d}\mathbf{r}
      + \sum_A Z_A \phi(\mathbf{R}_A) .

The first term enters the one-electron Hamiltonian as a matrix built by
``Ewald_Potential::IntegrateEwald``; the second is a scalar added to
:math:`E_0`. The minus sign on the electronic term is the electron charge,
and matches what ``AOMultipole::FillPotential`` does for the nuclei and for
external multipoles.

That sign is worth flagging because getting it wrong is nearly invisible.
For a neutral QM region the two terms differ only through the *shape* of the
density, not its total charge, so they very nearly cancel; flipping one
turns a near-cancellation into a near-doubling, and the result still looks
like a plausible small number. On a neutral methane in a methane
background, the difference is :math:`-13.5` meV against a correct
:math:`+0.13` meV.

The grid must be the one ``DFTEngine`` uses, since the potential values have
to land on exactly the quadrature points that form
:math:`\langle\rho|\phi\rangle`. ``PrepareEwaldPotentialGrid`` therefore
reads ``dftpackage.xtpdft.integration_grid``, not the ``grid_for_potential``
option — the latter governs the opposite direction, the QM density's
influence outward on the classical regions. The coupling between the two
grids is a constraint rather than a preference, and it costs nothing:
the Ewald contribution is converged at the ``medium`` default to below
:math:`10^{-10}` Ha.

Because the background is frozen, the QM geometry is fixed for the job, and
the QM density does not enter :math:`\phi`, the potential is evaluated once
per job and reused across inter-region SCF iterations.

The probe and Thole damping
===========================

:math:`\phi` at a point is obtained by placing a unit test charge there and
calling the same validated energy routines the classical channels use, so
that the result is in the same gauge by construction. That probe is a
``PolarSite``, whose constructor unavoidably assigns a polarizability from
the element table — which would otherwise feed a meaningless number into the
Thole damping of the induced-dipole term.

The field-point path therefore requests an **undamped** induced-dipole
interaction explicitly. This is not a preference: ``AOMultipole`` and
``DFTEngine::ExternalRepulsion`` already deliver induced dipoles to a QM
density undamped, so damping here would give one job two different
conventions for the same physical interaction depending on which route the
dipole arrived through.

Note that this cannot be expressed by giving the probe *no* polarizability.
The Thole factor would then be formed from :math:`u^3 = 0`, which the model
reads as complete overlap and damps maximally — the opposite of the intended
limit.

An analytic alternative, in which the background's moments are handed to the
DFT engine as operators rather than as a sampled potential, is partially
present in the code but not reachable. It is blocked on the
operator-centre derivatives a rank-1 source requires, which stock ``libint2``
builds do not provide.

.. _ewald-gauge:

The :math:`k = 0` term and the potential's gauge
================================================

The reciprocal sum omits the :math:`k = 0` term, which corresponds to
embedding the cell in a uniform neutralizing background. As a consequence
:math:`\phi` is defined only up to an additive constant :math:`\phi_0`, and
a region of net charge :math:`q` shifts by :math:`q\,\phi_0`.

For a neutral region this drops out identically:
:math:`\phi_0` multiplies :math:`Z_{\text{total}} - N_{\text{electrons}} = 0`
for the QM region, and :math:`\sum q = 0` for each neutral classical
segment. For a charged region it does not, and the resulting energy carries
the usual charged-periodic-cell convention. This is a physical statement
about the model, not a numerical artefact, and it is the same convention the
classical channel uses — every term in ``PotentialAt`` is evaluated through
the same energy routines, so the two halves of a QM/MM job share one gauge
by construction rather than by coincidence.


Practical notes
***************

.. _ewald-parameter-inheritance:

Parameter inheritance
=====================

**An** ``ewaldregion`` **takes every Ewald parameter from the background
checkpoint, and none from the job's own options.** :math:`\alpha`,
:math:`k_\text{max}`, :math:`r_\text{min}`, the field tolerance, the Thole
parameter, the screening factor, the shape and the box are all read from the
``ewald_parameters`` group written by ``ewaldbackground``. The job's XML
cannot override them, and does not offer the option.

This is deliberate. The background's induced dipoles were converged *with*
a particular splitting and particular cutoffs; the foreground's interaction
with that background is only meaningful if it is evaluated the same way. A
job that re-specified :math:`\alpha` would be adding a foreground computed
under one decomposition to a background computed under another, and the
:math:`\alpha` invariance of :ref:`ewald-alpha-invariant` would no longer
hold — silently, since each half would look internally consistent.

The practical consequence is that **changing any Ewald parameter means
re-running** ``ewaldbackground``. The checkpoint is the unit of
configuration, not the job file.

Choosing :math:`\alpha`
=======================

Since the total is independent of :math:`\alpha`, the choice is purely
one of cost. Real-space work scales as :math:`(s/\alpha)^3` with the
screening factor :math:`s`; reciprocal work scales as
:math:`k_\text{max}^3`, and the derived default :math:`k_\text{max} = 6\alpha`
makes that :math:`\alpha^3`. Balancing the two gives

.. math::

    \alpha_{\text{opt}} \approx \left(\frac{8\pi^3 N}{V^2}\right)^{1/6} ,

with :math:`N` the number of sites and :math:`V` the cell volume, which for
a few thousand sites in a few hundred thousand :math:`a_0^3` is of order
:math:`3\,\text{nm}^{-1}`.

If :math:`k_\text{max}` is set explicitly rather than derived, the balance
changes: reciprocal work becomes independent of :math:`\alpha` and larger
:math:`\alpha` is monotonically cheaper, up to the point where
:math:`k_\text{max}` is no longer large enough to converge the
:math:`\exp(-k^2/4\alpha^2)` weight. A ratio
:math:`k_\text{max}/\alpha = 6` corresponds to a truncation of
:math:`e^{-9}`; much beyond that is wasted work.

Cost
====

Evaluating the background potential over a DFT integration grid is the
dominant cost of a QM/MM job, and scales as the product of the number of
grid points, the number of background segments within the real-space cutoff,
and the number of :math:`k`-vectors. Those three numbers are logged before
the evaluation starts, along with periodic progress, so that a long run can
be distinguished from a hung one. The evaluation is parallelized over grid
points and happens once per job.

Diagnostics
===========

Two checks are worth running whenever this code is changed.

The :math:`\alpha` **scan** described in :ref:`ewald-alpha-invariant`: run
the same job against backgrounds converged at several :math:`\alpha`, with
all other convergence parameters held fixed, and confirm the total does not
move while the individual channels do.

The **with/without comparison**: run a job with and without the
``ewaldregion``, at otherwise identical settings, and difference the QM
energies at the first inter-region iteration, before any induced dipoles
exist. That difference is the QM–background interaction and nothing else. It
is the sharpest instrument available on this code path, and it is what
exposed the electron-charge sign described above — which the
zeroed-background test, the unit tests, and the :math:`\alpha` scan had all
missed.
