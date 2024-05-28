Theory and Background
#####################

Density-Functional Theory 
*************************

Excited state calculations require a reference ground state calculation within density-functional theory. VOTCA-XTP provides both an automated interface to the ORCA package [Neese:2012]_ and a lightweight internal DFT engine based on atom-centered Gaussian-type orbitals for method developing and testing. It solves the Kohn-Sham Equations for the molecular orbitals :math:`\phi_n^\textrm{KS}(\mathbf{r})` with orbital energies :math:`\varepsilon_n^\textrm{KS}$`

.. math::
    :label: equ:KS

    \left\{ -\frac{\hbar^2}{2m}\nabla^2 + V_\mathrm{ext}(\mathbf{r}) + V_\textrm{H}(\mathbf{r}) +V_\textrm{xc}(\mathbf{r})\right\}\phi_n^\textrm{KS}(\mathbf{r})  =\varepsilon_n^\textrm{KS} \phi_n^\textrm{KS}(\mathbf{r}) ,

where :math:`V_\textrm{ext}` is the external potential, :math:`V_\textrm{H}` the Hartree potential, and :math:`V_\textrm{xc}`` the exchange-correlation potential. VOTCA-XTP also contains functionality for projector-based-embedding DFT-in-DFT ground state calculations [Manby:2012]_, in which a chosen *active* subregion of a molecular system is embedded into an inactive one, reproducing the total energy of the full system ground state exactly.  

Many-Body Green's Functions and the Bethe-Salpeter Equation
***********************************************************

Using the ground-state reference, many-body Green's functions theory with the :math:`GW` approximation first calculates *single-particle excitations* (electron addition or removal) as solutions to the *quasiparticle equations*

.. math::
    :label: equ:QP

    \left\{ -\frac{\hbar^2}{2m}\nabla^2 + V_\textrm{ext}(\mathbf{r}) + V_\textrm{H}(\mathbf{r})\right\}\phi_n^\textrm{QP}(\mathbf{r}) + \int{\Sigma(\mathbf{r},\mathbf{r}',\varepsilon_n^\textrm{QP})\phi_n^\textrm{QP}(\mathbf{r}')d\mathbf{r}'} = \varepsilon_n^\textrm{QP} \phi_n^\textrm{QP}(\mathbf{r}) .

In place of the exchange-correlation potential in Eq.:eq:`equ:KS`, the energy-dependent self-energy operator :math:`\Sigma(\mathbf{r},\mathbf{r}',E) occurs in the QP equations. This operator is evaluated using the one-body Green's function in quasi-particle approximation

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



Hopping transport
*****************


