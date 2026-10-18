Polarizable Environments: Thole and Charge-Response Models
##########################################################

In a QM/MM job the classical environment of the QM region is divided into
regions (``staticregion``, ``polarregion``, ``ewaldregion``). Static regions
carry permanent multipoles only. Polar regions also respond to electric
fields, and this chapter describes how. A polar region can contain
molecules of two kinds, which can be mixed in one region:

* **Thole molecules**: inducible point dipoles on the atoms, coupled to each
  other inside and between molecules with Thole damping [Thole:1981]_. This
  has been XTP's model so far.

* **Charge-response molecules**: induced atomic charges from a
  charge-response matrix :math:`\Pi` of the molecule, plus residual induced
  dipoles, both smeared over Gaussians, with no couplings inside the
  molecule [Morita:1997]_ [DAvino:2014]_.

Both models give the same molecular polarizability tensor when they are
parametrized to it. They differ in *where* on the molecule the response
sits, which matters close to the QM region: at van der Waals contacts the
QM charge density probes the response atom by atom. In conjugated
molecules most of the in-plane polarizability is intramolecular charge flow
[Tsiper:2001]_, which induced charges describe directly and atomic dipoles
only effectively.

The same model, with the same smearing, is used in the ground-state QM/MM
loop and in the reaction field of environment-screened GW-BSE
(:doc:`embedded_gwbse`).


Common ground
*************

Sites and moments
=================

Every polar site carries the permanent multipoles of its ``.mps`` file (up
to quadrupoles) and a polarizability tensor :math:`\alpha_i` (the ``P``
line, in Å\ :sup:`3`, isotropic or as the six components
``xx xy xz yy yz zz``). The site responds with an induced dipole
:math:`\boldsymbol{\mu}_i`; a charge-response site in addition with an
induced charge :math:`\delta q_i`.

Fields and potentials
=====================

During an outer QM/MM iteration every other region adds its field (and, at
charge-response sites, its potential) to the sites of a polar region:

* permanent multipoles of static regions and of the other molecules of the
  same or other polar regions, as point multipoles;
* the induced moments of other polar regions;
* the QM region: the field of its electron density and nuclei;
* the periodic background of an ``ewaldregion`` (Thole molecules only; see
  `Limitations`_).

Two accumulators per site keep apart what is counted in the polar region's
energy (sources in outer regions) and what is counted in the other region's
energy (for example the QM region, whose Hamiltonian contains the induced
moments).

Energies
========

A polar region reports

.. math::

    E_\text{static} = E_\text{static-static} + E_\text{static,ext}, \qquad
    E_\text{polar} = E_\text{indu-indu} + E_\text{indu-stat}
                     + E_\text{internal} + E_\text{polar,ext},

with the interactions among induced moments of different sites
(:math:`E_\text{indu-indu}`), of induced moments with the permanent
multipoles of other molecules of the region (:math:`E_\text{indu-stat}`),
the internal (polarization) energy
:math:`E_\text{internal} = \tfrac12 \sum_i \boldsymbol{\mu}_i^T
\alpha_i^{-1} \boldsymbol{\mu}_i` (plus the charge term below for
charge-response molecules), and the interaction with the fields of outer
regions. At self-consistency :math:`E_\text{polar}` is the polarization
energy :math:`\tfrac12 x^T b` of the region in the fields it is counted
with.


The Thole model
***************

Each site responds to the total field, external plus that of all other
induced dipoles,

.. math::

    \boldsymbol{\mu}_i = \alpha_i \Big( \mathbf{E}_i^\text{ext}
       - \sum_{j \neq i} T_{ij}\, \boldsymbol{\mu}_j \Big)
    \quad\Leftrightarrow\quad
    (\alpha^{-1} + T)\, \boldsymbol{\mu} = \mathbf{E}^\text{ext},

where the sum runs over *all* other sites, inside the molecule as well. The
dipole tensor is damped at short range,

.. math::

    T_{ij} = \frac{\lambda_3}{r^3}\,\mathbb{1}
             - \frac{3\lambda_5}{r^5}\,\mathbf{r}\mathbf{r}^T, \qquad
    \lambda_3 = 1 - e^{-a u^3}, \quad
    \lambda_5 = 1 - (1 + a u^3)\, e^{-a u^3}, \quad
    a u^3 = \frac{a\, r^3}{\sqrt{\alpha_i \alpha_j}},

with :math:`a` the option ``exp_damp`` (default 0.39) and :math:`\alpha_i`
the largest eigenvalue of the site's polarizability tensor. Without the
damping, two point dipoles closer than :math:`(4 \alpha_i \alpha_j)^{1/6}`
polarize each other without bound (the polarization catastrophe).

Because the intramolecular couplings are part of the model, the site
polarizabilities are not the molecule's polarizability divided among its
atoms: ``xtp_tools -e molpol`` adjusts them until the coupled molecule
reproduces a target polarizability tensor, usually from a DFT calculation.

**Solver.** :math:`(\alpha^{-1} + T)\,\boldsymbol{\mu} = \mathbf{E}` is
solved by preconditioned conjugate gradients (``tolerance_dipole``,
``max_iter``), starting from the solution of the previous outer iteration;
a region of one molecule is solved directly.

**Coupling with the QM region.** The sites see the field of the QM
electrons integrated on the DFT grid at the site position and that of the
nuclei, both as points. The QM region sees the induced dipoles as point
dipoles in its Hamiltonian. In the reaction field of embedded GW-BSE the
sites are instead smeared over Gaussians of width
``site_width`` :math:`\cdot\,\alpha^{1/3}` to keep :math:`1 + R` positive
definite (:doc:`embedded_gwbse`).


The charge-response model
*************************

Response of one molecule
========================

A charge-response molecule :math:`I` responds to the external potential
:math:`\boldsymbol{\phi}_I` at its atoms with induced atomic charges and
to the field with residual induced dipoles,

.. math::

    \delta\mathbf{q}_I = -\Pi_I\, \boldsymbol{\phi}_I, \qquad
    \boldsymbol{\mu}_i = \alpha_i^\text{res}\, \mathbf{E}_i .

:math:`\Pi_I` (:math:`N \times N`, in e\ :sup:`2`/Hartree) is the
charge-response matrix of the molecule type: symmetric, positive
semidefinite, and charge conserving (each row sums to zero, so a uniform
potential induces nothing and every molecule stays neutral).
:math:`\alpha_i^\text{res}` carries what induced charges cannot: the
out-of-plane response of a planar molecule and a small in-plane residual.
In a uniform field :math:`\mathbf{F}` the potential is
:math:`\phi_k = -\mathbf{F}\cdot\mathbf{R}_k`, so the molecular
polarizability is

.. math::

    \alpha_\text{mol} = X^T \Pi\, X + \sum_i \alpha_i^\text{res},

with :math:`X` the atomic positions. The residual polarizabilities are fitted
so that this equals the DFT polarizability of the molecule exactly.

:math:`\Pi` is a linear-response property of the whole molecule and already
contains its own Coulomb screening. The induced charges and dipoles of one
molecule therefore do **not** interact with each other; they interact only
with other molecules and with the QM region.

Smearing
========

The induced charge and the residual dipole of a charge-response site are
Gaussians: the charge density :math:`\delta q\, g(\mathbf{r})` with
:math:`g(\mathbf{r}) = (\pi s^2)^{-3/2} e^{-r^2/s^2}`, and the dipole as
:math:`-\boldsymbol{\mu}\cdot\nabla g`. The width :math:`s` is the option
``cr_width`` (default 0.7 Å). Two Gaussians of widths :math:`s_i` and
:math:`s_j` interact through :math:`\operatorname{erf}(r/S)/r` with
:math:`S = \sqrt{s_i^2 + s_j^2}` and its derivatives,

.. math::

    B_0 = \frac{\operatorname{erf}(r/S)}{r}, \quad
    B_1 = \frac{B_0 - \frac{2}{\sqrt\pi\, S} e^{-r^2/S^2}}{r^2}, \quad
    B_2 = \frac{3 B_1 - \frac{4}{\sqrt\pi\, S^3} e^{-r^2/S^2}}{r^2},

so that, with :math:`\mathbf{r}` from the source to the receiving site, a
charge :math:`q` produces the potential :math:`q B_0` and the field
:math:`q B_1 \mathbf{r}`, and a dipole the potential
:math:`B_1\, \boldsymbol{\mu}\cdot\mathbf{r}` and the field
:math:`(B_2\, \mathbf{r}\mathbf{r}^T - B_1 \mathbb{1})\, \boldsymbol{\mu}`.
For :math:`S = 0` these are the bare Coulomb kernels. Thole sites are
points (:math:`s = 0`).

* Charge–charge and charge–dipole couplings always use these kernels.
* Dipole–dipole couplings between molecules, when at least one site is a
  charge-response site, use them as well (``cr_dipole_damping``
  ``gaussian``, the default) or the Thole damping above (``thole``).
  Between two Thole sites the Thole model is unchanged.
* Permanent multipoles act on charge-response sites as points. The
  Gaussian correction would be below :math:`3\cdot 10^{-4}` relative at the
  shortest intermolecular distances (1.7 Å for :math:`s = 0.7` Å).

The smearing is what makes the induced charges usable at contact
distances: a point charge right next to a hydrogen bond over-responds in the
same way a point dipole does. The QM/MM loop and the reaction field of
embedded GW-BSE use the same width, so the environment responds to the QM
density in the same way in both.

Linear system
=============

With the unknowns :math:`x = (\delta\mathbf{q}, \boldsymbol{\mu})` (induced
charges of all charge-response sites, induced dipoles of all sites), the
sources :math:`c = (-\boldsymbol{\phi}, \mathbf{E})` and
:math:`M = \operatorname{blockdiag}(\Pi_I, \alpha_i)`, the response of the
region is

.. math::

    x = M\,(c - T x)
    \quad\Leftrightarrow\quad
    (\mathbb{1} + M^{1/2}\, T\, M^{1/2})\, y = M^{1/2} c, \qquad
    x = M^{1/2} y,

where :math:`T` contains the couplings above: none inside a
charge-response molecule, the Thole couplings inside a Thole molecule, and
all couplings between molecules. :math:`T` is built so that :math:`(Tx)_q`
is the potential and :math:`(Tx)_\mu` minus the field of the induced moments.
The symmetric operator is positive definite for a stable environment and is
solved by conjugate gradients (``tolerance_dipole``, ``max_iter``). The null
space of :math:`\Pi` (charge conservation) is handled by :math:`M^{1/2}`.
A non-positive curvature during the iterations stops the job with a message
(polarization catastrophe). Regions without charge-response molecules use
the Thole solver unchanged.

The internal energy of the induced charges,
:math:`\tfrac12\, \delta\mathbf{q}^T \Pi^{+} \delta\mathbf{q}` with the
pseudo-inverse of :math:`\Pi`, is added to :math:`E_\text{internal}`.

Coupling with the QM region
===========================

The induced charge and dipole of a charge-response site enter the QM
Hamiltonian as the Gaussian charge and dipole above. The matrix elements are
three-centre Coulomb integrals of an AO pair with an s- and a p-type
primitive Gaussian at the site, normalized to a unit charge and a unit
dipole; their interaction with the QM nuclei uses the kernels :math:`B_0`
and :math:`B_1`. The potential and field of the QM density at the site,
seen by the same Gaussian, are traces of the density matrix with the same
integrals, and the nuclei again enter through :math:`B_0` and :math:`B_1`.
The coupling is therefore exactly symmetric, without grid errors. Only the
``xtp`` DFT package supports this; with other packages a job with
charge-response molecules stops at the first outer iteration.

Reaction field of embedded GW-BSE
=================================

With charge-response molecules among the explicit polar regions, the
reaction-field kernel over the auxiliary basis becomes

.. math::

    B = -(M^{1/2} C)^T\, (\mathbb{1} + M^{1/2} T M^{1/2})^{-1}\, (M^{1/2} C),
    \qquad C = (-\Phi^T, F^T),

with :math:`\Phi` the potential of each auxiliary function at the
charge-response sites (Gaussian-averaged with width ``cr_width``) and
:math:`F` its field at all sites (Thole sites smeared with ``site_width``
:math:`\cdot\,\alpha^{1/3}` as before). For Thole molecules alone it reduces
to :math:`-F (\alpha^{-1} + T)^{-1} F^T`. It is formed with a Cholesky
factor of the dense operator, as the Thole kernel. Charge-response molecules
cannot be part of a dielectric shell region (``shell_regions``), which
responds with :math:`\alpha/\varepsilon` per site.


Parametrization of charge-response molecules
********************************************

XTP does not compute :math:`\Pi`; it is the static linear response of the
molecule, obtained once per molecule type from a DFT calculation, like the
target polarizability of ``molpol``. The procedure the files described
below were validated with:

1. **Perturbations.** An external potential localized on atom :math:`k`
   with Hirshfeld weights [Hirshfeld:1977]_, :math:`v_k(\mathbf{r}) =
   w_k(\mathbf{r})` with :math:`\sum_k w_k = 1` (free spherical atoms in the
   same basis), plus three uniform fields.
2. **Response.** Analytic coupled-perturbed Kohn–Sham for all :math:`N + 3`
   perturbations at once; with diffuse functions in the basis (for example
   def2-SVPD or def2-TZVPD) and the functional of the QM region.
3. **Induced charges.** Measured with the same Hirshfeld weights:
   :math:`\Pi_{jk} = -\langle w_j | \chi | w_k \rangle`, which is symmetric,
   positive semidefinite and charge conserving by construction.
   (ESP-fitted charges, as in the original charge-response kernel
   [Morita:1997]_, give a non-symmetric, more ill-conditioned matrix with
   amplified short-wavelength modes; they were less accurate in the tests
   below.)
4. **Residual dipoles.**
   :math:`\alpha^\text{res} = \alpha_\text{DFT} - X^T \Pi X` from the field
   response, distributed over the atoms as one tensor scaled by the trace of
   the Hirshfeld atomic dipole polarizabilities. :math:`\alpha^\text{res}`
   must be positive definite.

Files
=====

* ``.crm``: the charge-response matrix in the atom order of the molecule
  type's ``.mps`` file. Lines starting with ``!`` are comments; then the
  number of atoms :math:`N`, :math:`N` lines ``element x y z`` (Å, the
  geometry of the fit, for information), and :math:`N` lines of :math:`N`
  numbers (e\ :sup:`2`/Hartree). On reading, XTP checks that it is
  symmetric, charge conserving and positive semidefinite to
  :math:`10^{-6}` relative, and removes the rounding of the file.
* ``.mps``: the permanent multipoles as for any polar molecule, with the
  residual polarizabilities :math:`\alpha_i^\text{res}` on the ``P`` lines,
  in the frame of the template geometry.
* map file: next to ``<multipoles_n>`` of the segment,

  .. code-block:: xml

      <multipoles_n>MP_FILES/water_cr_n.mps</multipoles_n>
      <chargeresponse_n>MP_FILES/water_cr_n.crm</chargeresponse_n>

  (and likewise for other states). The matrix is reordered to the order of
  ``mpoles`` in the fragments. Every atom of the molecule must be mapped: a
  part of :math:`\Pi` does not conserve charge. Segments mapped from an
  explicit ``.mps`` file name in the job file get no charge response.
* polar region options: ``cr_width`` and ``cr_dipole_damping`` (above);
  they apply to the charge-response molecules of the region.

Residual polarizabilities are in general anisotropic. Use ``map2md`` 0:
with ``map2md`` 1 the sites are placed on the MD atoms without rotating
site tensors.

Checkpoints store :math:`\Pi`, the induced charges and the options with the
region; checkpoints without them read as Thole molecules.


Accuracy and choice of model
****************************

The induction energy
:math:`P = \tfrac12 \langle \rho | \chi_\text{env} | \rho \rangle` of one
electron in a frontier orbital of a molecule A next to a molecule B (B the
environment), compared with B's DFT finite-field response (PBE0;
charge-response parameters at the same level):

.. list-table::
   :header-rows: 1

   * - A, B, contact
     - DFT (eV)
     - Thole
     - charge response, :math:`s = 0.7` Å
   * - pentacene pairs (bulk, 2.6–2.8 Å)
     - −0.10 to −0.17
     - +8 to +9 %
     - −4 to −10 %
   * - pyridine lone pair + H-bonded water (1.76 Å)
     - −0.318
     - +84 %
     - +5 %
   * - pyridine + 6 waters (coupled cluster of waters)
     - −0.449
     - +62 %
     - +5 %
   * - phenol + 6 waters
     - −0.178
     - +8 %
     - −3 %

For pentacene the Thole model was used with its ``molpol`` polarizabilities;
scaled to the DFT polarizability it is +15 to +17 %. Isotropic Thole sites
cannot follow the anisotropy of an acene (eigenvalues 21/41/68 Å\ :sup:`3`
vs. 18/37/91 Å\ :sup:`3` in DFT), while the charge-response model reproduces
the DFT tensor by construction. The width matters at hydrogen bonds (point
sites: +29 % for the pyridine pair) and hardly at all between π-stacked
acenes.

The charge-response model is the better choice for the molecules in direct
contact with the QM region when they are conjugated or hydrogen-bonded to
it. It costs one DFT linear-response calculation per molecule type; the
QM/MM loop itself is about as expensive as with Thole sites. Thole molecules
remain appropriate further out and can be mixed into the same region.


Limitations
***********

* Charge-response molecules need the ``xtp`` DFT package.
* They cannot be combined with an ``ewaldregion`` (its periodic sums act on
  point dipoles) or be part of a dielectric shell region of embedded
  GW-BSE.
* Whole molecules only (see `Files`_), and ``map2md`` 0 for anisotropic
  residual polarizabilities.
* The dense reaction-field kernel of embedded GW-BSE has the dimension
  :math:`N_\text{CR} + 3 N_\text{sites}`, slightly more than the
  :math:`3 N_\text{sites}` of the Thole model.
