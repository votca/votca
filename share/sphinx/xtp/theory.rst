Theory and Background
#####################

Density-Functional Theory 
*************************

Excited state calculations require a reference ground state calculation within density-functional theory. VOTCA-XTP provides both an automated interface to the ORCA package [Neese:2012]_ and a lightweight internal DFT engine based on atom-centered Gaussian-type orbitals for method developing and testing. It solves the Kohn-Sham Equations for the molecular orbitals :math:`\phi_n^\text{KS}(\mathbf{r})` with orbital energies :math:`\varepsilon_n^\textrm{KS}$`

.. math::
    :label: equ:KS

    \left\{ -\frac{\hbar^2}{2m}\nabla^2 + V_\mathrm{ext}(\mathbf{r}) + V_\textrm{H}(\mathbf{r}) +V_\textrm{xc}(\mathbf{r})\right\}\phi_n^\textrm{KS}(\mathbf{r})  =\varepsilon_n^\textrm{KS} \phi_n^\textrm{KS}(\mathbf{r}) ,

where :math:`V_\textrm{ext}` is the external potential, :math:`V_\textrm{H}` the Hartree potential, and :math:`V_\textrm{xc}` the exchange-correlation potential. VOTCA-XTP also contains functionality for projector-based-embedding DFT-in-DFT ground state calculations [Manby:2012]_, in which a chosen *active* subregion of a molecular system is embedded into an inactive one, reproducing the total energy of the full system ground state exactly.  

Many-Body Green's Functions and the Bethe-Salpeter Equation
***********************************************************

Using the ground-state reference, many-body Green's functions theory with the :math:`GW` approximation first calculates *single-particle excitations* (electron addition or removal) as solutions to the *quasiparticle equations*

.. math::
    :label: equ:QP

    \left\{ -\frac{\hbar^2}{2m}\nabla^2 + V_\textrm{ext}(\mathbf{r}) + V_\textrm{H}(\mathbf{r})\right\}\phi_n^\textrm{QP}(\mathbf{r}) + \int{\Sigma(\mathbf{r},\mathbf{r}',\varepsilon_n^\textrm{QP})\phi_n^\textrm{QP}(\mathbf{r}')d\mathbf{r}'} = \varepsilon_n^\textrm{QP} \phi_n^\textrm{QP}(\mathbf{r}) .

In place of the exchange-correlation potential in Eq.:eq:`equ:KS`, the energy-dependent self-energy operator :math:`\Sigma(\mathbf{r},\mathbf{r}',E)` occurs in the QP equations. This operator is evaluated using the one-body Green's function in quasi-particle approximation

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


As the correction :math:`\Delta \varepsilon_n^{GW}` itself depends on :math:`\varepsilon_n^\textrm{QP}`, Eq.:eq:`equ:theory:gw_sc` needs to be solved self-consistently.

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
Here, :math:`K^\text{x}` is the repulsive exchange interaction originating from the bare Coulomb term :math:`v_\mathrm{C}`, while the direct interaction :math:`K^\text{d}` contains the attractive, but screened, interaction :math:`W` between electron and hole, causing the binding of the electron-hole pair. In Eq.:eq:`equ:theory:Kd` it is assumed that the dynamic properties of :math:`W(\omega)` are negligible, and the computationally less demanding static approximation :math:`\omega=0` is employed. 




Hopping transport
*****************


