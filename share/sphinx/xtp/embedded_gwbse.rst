Embedded GW-BSE: Screening by a Polarizable Environment
#######################################################

An electronic excitation of a molecule in a solid is screened by its
surroundings. When the excitation rearranges charge, the neighbouring
molecules polarize in response, and that response acts back on the
excitation. For charged excitations (electron or hole) the effect is of the
order of an electronvolt; for neutral excitations it is smaller, but large
enough to reorder states or to shift them by more than the accuracy one
expects of GW-BSE.

VOTCA-XTP offers two ways to include this response in a QM/MM job with
polarizable regions:

* **Iteratively**, by treating the excited state like a ground state: the
  QM/MM loop converges the environment's induced dipoles self-consistently
  against the density of the excited state, with a GW-BSE calculation in
  every outer iteration. This is the default for a job whose QM region is
  in an excited state.

* **Embedded**, by folding the environment's response into the screened
  Coulomb interaction of GW-BSE itself, following
  [Li:2016]_ and [Li:2018]_. The QM/MM loop converges only the ground state;
  one GW-BSE calculation then runs with the environment as an additional,
  static screening medium. This is switched on with
  ``environment_screening`` in the ``qmregion`` options.

The embedded scheme is cheaper, since it needs one GW-BSE calculation
instead of one per outer iteration. It also contains two contributions that
the iterative scheme cannot represent (see :ref:`embedded-vs-iterative`).
The ``radial_dielectric`` correction of the legacy code is available within
it as shell regions.


The reaction field
******************

The environment is the set of polarizable sites of the job's polar
regions: Thole-damped inducible point dipoles [Thole:1981]_. A change
:math:`\delta\rho` of the QM charge density produces a field
:math:`\mathbf{E} = F\,\delta\rho` at the sites, the sites respond with
induced dipoles :math:`\boldsymbol{\mu} = A^{-1}\mathbf{E}`, where

.. math::

    A = \alpha^{-1} + T_\text{Thole}

is the same operator the polar region solves in the QM/MM loop, and the
induced dipoles act back on the QM subsystem. The resulting additional
interaction between two QM charge distributions is the reaction field

.. math::
    :label: equ:embedded:vreac

    v_\text{reac} = -F A^{-1} F^{T},

and the QM electrons interact through

.. math::

    u = v + v_\text{reac}

instead of the bare Coulomb interaction :math:`v`. :math:`v_\text{reac}` is
negative semidefinite: the environment screens.

In the resolution-of-identity representation used by GW-BSE, a product
density :math:`\rho_{mn}` is expanded in auxiliary functions,
:math:`F` becomes the field at the sites of each auxiliary function taken as
a charge density, and :eq:`equ:embedded:vreac` turns into a matrix
:math:`B = -F A^{-1} F^{T}` over the auxiliary basis. In the metric of the
stored three-centre integrals :math:`M_{mn}` (in which the bare :math:`v` is
the identity) it becomes

.. math::

    R = T^{T} B\, T, \qquad
    (mn|v_\text{reac}|kl) = M_{mn}\, R\, M_{kl}^{T},

with :math:`T` the inverse square root of the auxiliary Coulomb metric that
is folded into :math:`M`. :math:`R` is built once per job: it depends only
on the QM atoms, their auxiliary basis and the polar sites, not on the
orbitals.

:math:`A` is assembled and factorized densely, which is far cheaper than
iterating a linear solve for every auxiliary function, and :math:`B` is
formed from its Cholesky factor, so it is symmetric and negative
semidefinite by construction. The field :math:`F` is obtained from
nuclear-attraction (or, for smeared sites, two-centre Coulomb) integrals
of each auxiliary function by a four-point central difference.


Quasiparticle energies (GW)
***************************

With the environment, the screened interaction is

.. math::

    W = \left[u^{-1} - \chi_0\right]^{-1},

which is split as :math:`W = v + v_\text{reac} + (W - u)`:

* :math:`v` gives the exchange self-energy :math:`\Sigma^\text{x}`, as
  without an environment.

* :math:`v_\text{reac}` gives a static COH+SEX self-energy,

  .. math::

      \Sigma^\text{reac}_{mn} = \tfrac{1}{2} \sum_{l} s_l\,
          (ml|v_\text{reac}|ln), \qquad
      s_l = \begin{cases} -1 & l \text{ occupied} \\
                          +1 & l \text{ virtual,} \end{cases}

  that is, a Coulomb hole :math:`\tfrac{1}{2}\sum_l (ml|v_\text{reac}|ln)`
  plus a screened exchange :math:`-\sum_{l\in\text{occ}} (ml|v_\text{reac}|ln)`.
  It is part of the correlation (it belongs to :math:`W - v`) and enters
  wherever :math:`\Sigma^\text{c}` does: the quasiparticle equation, evGW
  and QSGW. It is printed as ``S-R`` next to ``S-X`` and ``S-C``, and the
  log lists its diagonal together with the resulting HOMO, LUMO and gap
  shifts.

* :math:`W - u` is the correlation of the QM electrons, which now interact
  through :math:`u` and are screened by the environment as well. With
  :math:`S = (1 + R)^{1/2}`,

  .. math::

      W = S \left[1 - S\chi_0 S\right]^{-1} S ,

  so the unchanged RPA and correlation machinery runs on the *dressed*
  integrals :math:`M S`.

For an electron added to or removed from the molecule,
:math:`\Sigma^\text{reac}` contains the polarization energy of the
environment by that charge. The HOMO and LUMO quasiparticle energies
therefore include it directly, which is how the embedded scheme provides
hole and electron site energies (see :ref:`embedded-usage`).


Neutral excitations (BSE)
*************************

The BSE kernel becomes

.. math::

    K = K^\text{x}[v] + K^\text{reac}[v_\text{reac}] + K^\text{d}[W_\text{tot}],
    \qquad
    W_\text{tot} = \left[(1 + R)^{-1} + \varepsilon - 1\right]^{-1},

with :math:`\varepsilon` the RPA dielectric matrix of the QM subsystem in
the same metric.

* :math:`K^\text{d}[W_\text{tot}]` is the direct (electron-hole) term with
  the environment-screened interaction. It describes the environment's
  response to the charge that the excitation rearranges: a state-specific
  effect, in the language of solvation models [Duchemin:2018]_.

* :math:`K^\text{reac} = (ia|v_\text{reac}|jb)` has the index structure of
  the exchange term. It is the environment's response to the *transition*
  density, i.e. a linear-response effect, and it mainly moves bright
  singlets. Triplets have no exchange term and therefore no
  :math:`K^\text{reac}`.

With ``include_kreac`` (the default), the whole kernel runs on dressed
integrals, where :math:`K^\text{x} + K^\text{reac}` and
:math:`K^\text{d}[W_\text{tot}]` come out of the unchanged BSE operators.
Without it, the integrals stay bare and :math:`W_\text{tot}` is
diagonalized in place of :math:`\varepsilon`, so that :math:`K^\text{x}`
keeps the bare :math:`v`. In both cases the first-order shift
:math:`\langle K^\text{reac}\rangle = 2\,(X+Y)^{T} K^\text{reac} (X+Y)` of
each singlet is logged after the BSE is solved, as an estimate of the
linear-response contribution. With ``include_kreac`` the exchange column of
the state analysis is labelled ``<K_x+K_reac>``.


Site widths and stability
*************************

For a charge density and inducible points *outside* it, the reaction
energy :math:`\langle\rho|v_\text{reac}|\rho\rangle` is bounded by the
dielectric limit, above :math:`-\langle\rho|v|\rho\rangle`. In the metric of
:math:`M` this means that all eigenvalues of :math:`R` lie above
:math:`-1`, i.e. that :math:`1 + R` is positive definite. A point dipole
*inside* the tail of a diffuse auxiliary function is not bounded in this
way, and in dense molecular solids the closest sites of neighbouring
molecules do sit in those tails. An eigenvalue of :math:`R` at or below
:math:`-1` would mean that the environment screens a charge fluctuation by
more than the fluctuation itself; the dressing :math:`S = (1+R)^{1/2}` then
does not exist.

Each site therefore responds to the QM field averaged over a normalized
Gaussian of width

.. math::

    R_j = w\,\alpha_{\text{iso},j}^{1/3},
    \qquad \alpha_{\text{iso},j} = \tfrac{1}{3}\operatorname{tr}\alpha_j ,

with :math:`w` the option ``site_width`` -- the length scale that Thole
damping itself uses. Outside the auxiliary functions a spherical charge
acts as a point, so this changes the response only where point sites
over-respond. In a production QM/MM job (55 QM atoms, def2-TZVP with
aux-def2-TZVP, 3900 Thole sites) the lowest eigenvalue of :math:`R` was
:math:`-1.056` with point sites, carried by three carbon atoms of a
neighbouring C\ :sub:`60` at 3.05--3.3 Å; with the default
``site_width`` of 0.5 it was :math:`-0.755`, while the reaction energies of
the HOMO, LUMO and HOMO-LUMO densities changed by :math:`10^{-4}` relative.

Because :math:`R` does not depend on the orbitals, a job checks
:math:`1 + R` *before* the QM/MM loop starts. If it is not positive
definite, the job fails right away with a report of the lowest modes: where
their charge sits on the QM side and which polar sites carry their
reaction, with distances.


Shell regions
*************

Polar regions listed in ``shell_regions`` respond without mutual
induction, each site with the polarizability
:math:`\alpha_j/\varepsilon_\text{shell}`:

.. math::

    B_\text{shell} = -\sum_{j} F_j\,
        \frac{\alpha_j}{\varepsilon_\text{shell}}\, F_j^{T} .

This needs no linear solve, so it is a cheap treatment of an outer shell
beyond the explicit Thole region. The factor
:math:`1/\varepsilon_\text{shell}` (``shell_dielectric``) stands in for the
mutual induction that is left out. It is the ``radial_dielectric``
correction of the legacy code in this framework. Because
the sum runs over actual sites, the geometry is that of the morphology (a
slab simply has no sites in the vacuum), which is why this rather than an
analytic continuum term is used for the tail. The shell and the explicit
regions do not couple to each other.

The explicit (non-shell) polar regions respond together, as one coupled
Thole system, and must therefore use the same ``exp_damp``.

An ``ewaldregion`` contributes the electrostatic potential of the periodic
background to the ground state, as in any QM/MM job (see
:doc:`ewald_embedding`). It is frozen and does not respond to the
excitation, so it is not part of :math:`R`.


.. _embedded-vs-iterative:

Relation to the iterative scheme
********************************

With a linear environment response, the total energy of the QM/MM system
in a state :math:`X` with QM density :math:`\rho_X` is, at
self-consistency,

.. math::

    E_X = E_\text{QM}[\rho_X] + \langle\rho_X|\phi_\text{p}\rangle
        + \tfrac{1}{2}\langle\rho_X|v_\text{reac}|\rho_X\rangle ,

with :math:`\phi_\text{p}` the potential of the permanent multipoles. The
excitation energy of the iterative scheme is then

.. math::

    \Delta E = E_X - E_g = \Omega_g
        + \tfrac{1}{2}\langle\Delta\rho|v_\text{reac}|\Delta\rho\rangle ,

where :math:`\Omega_g` is the excitation energy in the environment
polarized by the ground state and :math:`\Delta\rho = \rho_X - \rho_g` is
the density change of the excitation. The iterative scheme thus sees the
environment's response to the *mean* density change of the excitation.

The embedded scheme starts from the same ground state, so its
:math:`\Omega_g` is the same, and through :math:`K^\text{d}[W_\text{tot}]`
it contains the mean-field response as well. In addition, it contains

* the response to the correlated motion of electron and hole: the
  environment responds to where the electron and the hole are relative to
  each other, not only to their averaged densities; and

* :math:`K^\text{reac}`, the response to the transition density.

Neither can be represented by a density that is iterated against the
environment. Conversely, both schemes treat the environment's response as
static (instantaneous); see [Amblard:2024]_ for when a dynamically
polarizable environment matters.

The same distinction applies to charged excitations. In a charged QM/MM
job (state ``e`` or ``h``) the environment polarizes against the averaged
density of the extra electron or hole. The static self-energy
:math:`\Sigma^\text{reac}` responds to the carrier as an instantaneous
point charge in each orbital it occupies. The two can differ noticeably
for orbitals spread over several separated parts of a molecule.


.. _embedded-usage:

Usage
*****

The embedded scheme is switched on per QM region:

.. code-block:: xml

    <qmregion>
      <id>0</id>
      <state>jobfile</state>
      <segments>jobfile</segments>
      <dftpackage> ... </dftpackage>
      <gwbse> ... </gwbse>
      <statetracker> ... </statetracker>
      <environment_screening>
        <include_kreac>true</include_kreac>
        <shell_regions>2</shell_regions>
        <shell_dielectric>4.0</shell_dielectric>
        <site_width>0.5</site_width>
      </environment_screening>
    </qmregion>

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - option
     - default
     - meaning
   * - ``include_kreac``
     - ``true``
     - Include :math:`K^\text{reac}` in the BSE kernel. Without it,
       :math:`K^\text{reac}` is only estimated to first order and logged.
   * - ``shell_regions``
     - (none)
     - Ids of polar regions that respond as an uncoupled shell, with
       polarizabilities :math:`\alpha/\varepsilon_\text{shell}`.
   * - ``shell_dielectric``
     - ``4.0``
     - :math:`\varepsilon_\text{shell}` for the shell regions.
   * - ``site_width``
     - ``0.5``
     - Width of the polar sites as seen by the QM density, in units of
       :math:`\alpha_\text{iso}^{1/3}`; 0 gives point sites.

The state of the QM region must be a GW-BSE state:

* an exciton (``s1``, ``t2``, ...), tracked with the ``statetracker`` as
  usual, or
* a quasiparticle level (``pqpN`` or ``dqpN``, with ``N`` the HOMO or
  LUMO index), for hole and electron site energies.

Charged DFT states (``e``, ``h``) are refused. Their extra charge would be
polarized in the QM/MM loop and screened a second time in GW. Unrestricted
excitons and truncated (DFT-in-DFT) active regions are not supported.

A job proceeds as follows:

#. The screening environment is assembled from the job's polar regions,
   and :math:`1 + R` is checked.
#. The QM/MM loop converges the **ground state**: DFT and the polar
   regions, as in an ``n`` job.
#. One GW-BSE calculation runs on top of it, screened by the environment.
   The tracked state's energy replaces the QM region's energy.
#. The job writes ``checkpoint_screened.hdf5`` and reports, in its
   output, the ground-state total energy ``E_ground`` and the site energies
   ``screened_site_energies``, all relative to its own ground state, in eV:

   * ``h`` :math:`= -\varepsilon_\text{HOMO}` and
     ``e`` :math:`= \varepsilon_\text{LUMO}`, the quasiparticle energies
     (diagonalized ones if available),
   * ``s`` or ``t``, the excitation energy of the tracked exciton, if the
     job's state is one.

A screened job therefore needs no separate ground-state job. When writing
jobs (``-j write``) no ``n`` job is added for the segments, and ``-j read``
takes the site energies of screened jobs directly. Several screened jobs of
the same segment (an ``s1`` and a quasiparticle job, say) report the same
``e`` and ``h``; ``-j read`` checks that they agree to 1 meV.

The auxiliary basis of the screening is the one GW-BSE uses: the DFT
auxiliary basis if one is set, otherwise ``gwbse.auxbasisset``, otherwise
``aux-<basisset>``.


Example
*******

For a dicyanovinyl-oligothiophene donor in a morphology with C\ :sub:`60`
(55 QM atoms,
def2-TZVP/aux-def2-TZVP, evGW, full BSE; polar region of 2.4 nm with 68
segments and 3900 Thole sites, ``exp_damp`` 0.39; optionally a periodic
Ewald background of 1212 segments; ``site_width`` 0.5, no shell regions),
the lowest singlet excitation energy (eV) is

.. list-table::
   :header-rows: 1
   :widths: 46 22 32

   * - scheme
     - cut-off
     - with Ewald background
   * - QM only (vacuum)
     - 2.378
     -
   * - iterative QM/MM, :math:`E_{s_1} - E_n`
     - 2.319
     - 2.286
   * - embedded, without :math:`K^\text{reac}`
     - 2.218
     - 2.186
   * - embedded, with :math:`K^\text{reac}`
     - 2.102
     - 2.068

The Ewald background shifts all three schemes by the same
:math:`-0.033` eV: it acts through the ground state only. Relative to the
QM-only value (cut-off environment), the ground-state embedding contributes
:math:`+0.027` eV and the mean-field response to the density change
:math:`-0.087` eV in the iterative and :math:`-0.083` eV in the embedded
scheme. These two are common to both schemes. The embedded scheme adds :math:`-0.104` eV from the correlated electron-hole pair and
:math:`-0.116` eV from :math:`K^\text{reac}`, which the logged first-order
estimate reproduces to about 10%.
