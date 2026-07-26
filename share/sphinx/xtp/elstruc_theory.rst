Electronic Structure Theories
#############################

Density-Functional Theory 
*************************

Excited state calculations require a reference ground state calculation within density-functional theory. VOTCA-XTP provides both an automated interface to the ORCA package [Neese:2012]_ and a lightweight internal DFT engine based on atom-centered Gaussian-type orbitals for method developing and testing. It solves the Kohn-Sham Equations for the molecular orbitals :math:`\phi_n^\textrm{KS}(\mathbf{r})` with orbital energies :math:`\varepsilon_n^\textrm{KS}`

.. math::
    :label: equ:KS

    \left\{ -\frac{\hbar^2}{2m}\nabla^2 + V_\mathrm{ext}(\mathbf{r}) + V_\textrm{H}(\mathbf{r}) +V_\textrm{xc}(\mathbf{r})\right\}\phi_n^\textrm{KS}(\mathbf{r})  =\varepsilon_n^\textrm{KS} \phi_n^\textrm{KS}(\mathbf{r}) ,

where :math:`V_\textrm{ext}` is the external potential, :math:`V_\textrm{H}` the Hartree potential, and :math:`V_\textrm{xc}` the exchange-correlation potential. VOTCA-XTP also contains functionality for projector-based-embedding DFT-in-DFT ground state calculations [Manby:2012]_, in which a chosen *active* subregion of a molecular system is embedded into an inactive one, reproducing the total energy of the full system ground state exactly.

For closed-shell systems, :eq:`equ:KS` is solved once for a common set of spatial orbitals doubly occupied by electrons of opposite spin (restricted Kohn-Sham, RKS). Open-shell systems -- radicals, cations, anions, or any system with an odd number of electrons or a triplet or higher-multiplicity ground state -- instead require separate equations for the :math:`\alpha` and :math:`\beta` spin channels, each with its own set of orbitals, orbital energies, and occupied-orbital count (unrestricted Kohn-Sham, UKS). The internal DFT engine automatically switches to the UKS formalism whenever the numbers of :math:`\alpha` and :math:`\beta` electrons differ, as is always the case, for example, for a molecule with an odd total charge or an unpaired electron.

Self-consistent field convergence
==================================

Because :math:`V_\textrm{H}` and :math:`V_\textrm{xc}` in :eq:`equ:KS` themselves depend on the orbitals being solved for, the Kohn-Sham equations must be solved iteratively: a trial density produces a Fock (Kohn-Sham) matrix, which is diagonalized to give new orbitals and a new density, and this cycle repeats until the change in energy and density falls below a set of convergence thresholds. Convergence of this self-consistent field (SCF) procedure is normally accelerated using the direct inversion of the iterative subspace (DIIS) [Pulay:1980]_ and its energy-based variant (A-DIIS) [Hu:2010]_, which extrapolate a new Fock matrix guess from a short history of previous iterations rather than using the current iteration's own Fock matrix directly. In addition, level shifting of the virtual orbitals and linear mixing of the density with that of the previous iteration are available to further stabilize convergence.

Open-shell systems with a small HOMO-LUMO gap, or in which the physically correct electronic state involves a specific localization of charge or spin density that is not obviously the lowest-energy solution to plain SCF, are systematically more difficult to converge with DIIS-based methods alone -- a well-documented limitation of DIIS-family extrapolation rather than a peculiarity of any one implementation. For such cases, the internal DFT engine falls back to direct energy minimization in orbital-rotation space once DIIS-based convergence has repeatedly failed, or once the SCF error is observed to be improving too slowly on average over a trailing window of iterations. This fallback solves an augmented-Hessian eigenvalue problem [Helmich-Paris:2021]_ -- the same family of algorithm used, for example, in ORCA's Trust-Region Augmented-Hessian (TRAH) SCF converger -- within a trust-region framework that only accepts a step once it is confirmed, after the fact, to actually lower the energy by an amount consistent with the quadratic model used to propose it. Further implementation details are given in :doc:`implementation_notes`.

Constrained Density-Functional Theory
**************************************

Ordinary Kohn-Sham DFT finds the ground-state density by minimizing the total energy functional :math:`E[n]` over all densities integrating to the correct total number of electrons. For many systems of interest -- localized charge or spin states relevant to electron- and hole-transfer reactions, diabatic states for computing electronic couplings, or simply a well-defined starting point for systems with a near-degenerate manifold of low-lying electronic configurations -- the physically relevant state is not necessarily the unconstrained minimum. Self-interaction error in approximate exchange-correlation functionals in particular tends to systematically favor spurious delocalization of charge or spin density over two or more fragments, even when the physically correct picture is a localized one.

Constrained density-functional theory (CDFT) [Wu:2005]_ addresses this by minimizing :math:`E[n]` subject to an additional constraint: that some population :math:`N_c`, obtained by integrating the density against a chosen weight function :math:`w_c(\mathbf{r})`, take on a specified target value,

.. math::
    :label: equ:cdft_constraint

    \int w_c(\mathbf{r})\, n(\mathbf{r})\, d\mathbf{r} = N_c .

The weight function :math:`w_c(\mathbf{r})` defines which part of the system the constraint applies to. VOTCA-XTP uses Hirshfeld (stockholder) partitioning [Hirshfeld:1977]_, in which each atom :math:`A` is assigned a weight

.. math::
    :label: equ:hirshfeld_weight

    w_A(\mathbf{r}) = \frac{n_A^{\textrm{at}}(\mathbf{r})}{\sum_B n_B^{\textrm{at}}(\mathbf{r})}

built from a superposition of free-atom reference densities :math:`n_A^{\textrm{at}}`, and the constraint weight for a group (fragment) of atoms is the sum of the weights of its constituent atoms. This has been found to give a more physically faithful charge partitioning than simpler, purely geometric schemes for the class of charge- and spin-localization problems CDFT is typically used for [Ahart:2022]_.

The constraint in :eq:`equ:cdft_constraint` is enforced via a Lagrange multiplier :math:`V_c`, leading to the constrained functional

.. math::
    :label: equ:cdft_functional

    W[n, V_c] = E[n] + V_c\left(\int w_c(\mathbf{r})\, n(\mathbf{r})\, d\mathbf{r} - N_c\right) .

The constrained ground state is the density that minimizes :math:`W` for the multiplier value :math:`V_c` that itself maximizes :math:`W` -- a saddle-point problem, though one with a unique solution [Wu:2005]_. In practice, this is solved as two nested loops: an inner Kohn-Sham SCF loop for fixed :math:`V_c`, in which the constraint enters as an additional, weight-function-dependent potential added directly to the Kohn-Sham Hamiltonian in :eq:`equ:KS`,

.. math::
    :label: equ:cdft_potential

    V_\textrm{c}(\mathbf{r}) = V_c\, w_c(\mathbf{r}) ,

and an outer bisection loop over :math:`V_c` itself, adjusting the multiplier until the population in :eq:`equ:cdft_constraint` matches the requested target. Since the inner SCF for a nearby value of :math:`V_c` typically provides a good starting guess for the next, the outer loop warm-starts each new inner SCF from the orbitals of the previous one.

Because CDFT constrains a population defined via a weight function that depends on the nuclear positions (both directly, through the definition of the atomic fragments, and indirectly, through the free-atom reference densities entering :eq:`equ:hirshfeld_weight`), the force acting on the nuclei is not given by the ordinary Hellmann-Feynman expression alone. Differentiating the constrained functional with respect to a nuclear coordinate yields additional terms describing how the weight function itself, the associated real-space integration grid, and the basis functions centered on the moving atoms all change with the nuclear positions [Oberhofer:2009]_. VOTCA-XTP evaluates these terms explicitly, so that CDFT forces are available for geometry optimization and molecular dynamics on the constrained potential energy surface, not just single-point energies.

Many-Body Green's Functions and the Bethe-Salpeter Equation
***********************************************************

Using the ground-state reference, many-body Green's functions theory with the :math:`GW` approximation first calculates *single-particle excitations* (electron addition or removal) as solutions to the *quasiparticle equations*

.. math::
    :label: equ:QP

    \left\{ -\frac{\hbar^2}{2m}\nabla^2 + V_\textrm{ext}(\mathbf{r}) + V_\textrm{H}(\mathbf{r})\right\}\phi_n^\textrm{QP}(\mathbf{r}) + \int{\Sigma(\mathbf{r},\mathbf{r}',\varepsilon_n^\textrm{QP})\phi_n^\textrm{QP}(\mathbf{r}')d\mathbf{r}'} = \varepsilon_n^\textrm{QP} \phi_n^\textrm{QP}(\mathbf{r}) .

In place of the exchange-correlation potential in :eq:`equ:KS`, the energy-dependent self-energy operator :math:`\Sigma(\mathbf{r},\mathbf{r}',E)` occurs in the QP equations. This operator is evaluated using the one-body Green's function in quasi-particle approximation

.. math::
    :label: equ:Green

    G(\mathbf{r},\mathbf{r}',\omega) = \sum_n{\frac{\phi_n(\mathbf{r})\phi_n^*(\mathbf{r}')}{\omega-\varepsilon_n+i0^+\textrm{sgn}(\varepsilon_n -\mu)}}

as

.. math::
    :label: equ:sigma
    
    \Sigma(\mathbf{r},\mathbf{r}',E) = \frac{i}{2\pi} \int{e^{-i\omega 0^+}G(\mathbf{r},\mathbf{r}',E-\omega)W(\mathbf{r},\mathbf{r}',\omega)\,d\omega},

where :math:`W` denotes the dynamically screened Coulomb interaction. Assuming that :math:`\phi^\textrm{QP}_n\approx \phi^\textrm{KS}_n`, the quasiparticle energies can be evaluated perturbatively according to

.. math::
    :label: equ:theory:gw_sc

    \varepsilon_n^\textrm{QP}= \varepsilon_n^\textrm{KS} + \Delta \varepsilon_n^{GW} = \varepsilon_n^\textrm{KS} + \left\langle\phi^\textrm{KS}_n\left\vert \Sigma(\varepsilon_n^\textrm{QP})-V_\text{xc} \right\vert\phi^\textrm{KS}_n\right\rangle .

As the correction :math:`\Delta \varepsilon_n^{GW}` itself depends on :math:`\varepsilon_n^\textrm{QP}`, :eq:`equ:theory:gw_sc` needs to be solved self-consistently.

Neutral excitations with a conserved number of electrons can be obtained from the Bethe-Salpeter Equation (BSE) by expressing coupled electron-hole amplitudes of excitation :math:`S` in a product basis of single-particle orbitals, i.e., 

.. math:: 
    :label: equ:bsewf
    
    \chi_S(\mathbf{r}_\textrm{e},\mathbf{r}_\textrm{h})=\sum_{v}^{\mathrm{occ}}\sum_c^{\mathrm{unocc}}A_{vc}^S\phi_{c}(\mathbf{r}_\textrm{e})\phi^*_{v}(\mathbf{r}_\textrm{h})+B_{vc}^S\phi_{v}(\mathbf{r}_\textrm{e})\phi^{*}_{c}(\mathbf{r}_\textrm{h}),
  
where :math:`\mathbf{r}_\textrm{e}` (:math:`\mathbf{r}_\textrm{h}`) is for the electron (hole) coordinate and :math:`A_{vc}` (:math:`B_{vc}`) are the expansion coefficients of the excited state wave function in terms of resonant (anti-resonant) transitions between occupied :math:`v` and unoccupied :math:`c` states, respectively. In this basis, the BSE turns into an effective two-particle Hamiltonian problem of the form

.. math::
    :label: equ:theory:bseeigenvalue

    \begin{pmatrix}
                                    \underline{\mathbf{H}}^{\text{res}}&\underline{\mathbf{K}} \\
                                    -\underline{\mathbf{K}} & -\underline{\mathbf{H}}^{\text{res}}
                                   \end{pmatrix}
                                   \begin{pmatrix}
                                   \mathbf{A}^S\\ \mathbf{B}^S
                                   \end{pmatrix}
                                   =\Omega_S
                                   \begin{pmatrix}
                                   \mathbf{A}^S\\ \mathbf{B}^S 
                                   \end{pmatrix}.

Specifically, the matrix elements of the blocks :math:`\underline{\mathbf{H}}^{\text{res}}` and :math:`\underline{\mathbf{K}}` are calculated as

.. math::
    :label: equ:BSEelements

    \begin{align}
    H^{\text{res}}_{vc,v'c'}&=D_{vc,v'c'}+\eta K^\mathrm{x}_{vc,v'c'}+K^\mathrm{d}_{vc,v'c'}\\
     K_{cv,v'c'}&=\eta K^\mathrm{x}_{cv,v'c'}+K^\mathrm{d}_{cv,v'c'}\, ,
    \end{align}

with 

.. math::
    :label: equ:theory:Kd

    \begin{align}
        D_{vc,v'c'}&=(\varepsilon_c-\varepsilon_v)\delta_{vv'}\delta_{cc'},\\
        K^\text{x}_{vc,v'c'}&=\iint  \phi_c^*(\mathbf{r}_\textrm{e})\phi_v(\mathbf{r}_\textrm{e})v_{\mathrm{C}}(\mathbf{r}_\textrm{e},\mathbf{r}_\textrm{h}) \phi_{c'}(\mathbf{r}_\textrm{h})\phi_{v'}^*(\mathbf{r}_\textrm{h}) d^3\mathbf{r}_\textrm{e} d^3\mathbf{r}_\textrm{h}\\
        K^\text{d}_{vc,v'c'}&=-\iint
               \phi_c^*(\mathbf{r}_\textrm{e})\phi_{c'}(\mathbf{r}_\textrm{e})W(\mathbf{r}_\textrm{e},\mathbf{r}_\textrm{h},\omega=0) \phi_v(\mathbf{r}_\textrm{h})\phi_{v'}^*(\mathbf{r}_\textrm{h})d^3\mathbf{r}_\textrm{e} d^3\mathbf{r}_\textrm{h}
               \, .
    \end{align}

and :math:`\eta=2` (:math:`\eta=0`) for singlet (triplet) excitations.
Here, :math:`K^\text{x}` is the repulsive exchange interaction originating from the bare Coulomb term :math:`v_\mathrm{C}`, while the direct interaction :math:`K^\text{d}` contains the attractive, but screened, interaction :math:`W` between electron and hole, causing the binding of the electron-hole pair. In :eq:`equ:theory:Kd` it is assumed that the dynamic properties of :math:`W(\omega)` are negligible, and the computationally less demanding static approximation :math:`\omega=0` is employed. 


