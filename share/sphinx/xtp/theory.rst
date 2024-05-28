Theory and Background
#####################

Density-Functional Theory 
*************************

Excited state calculations require a reference ground state calculation within density-functional theory. VOTCA-XTP provides both an automated interface to the ORCA package [Neese:2012]_ and a lightweight internal DFT engine based on atom-centered Gaussian-type orbitals for method developing and testing. It solves the Kohn-Sham Equations for the molecular orbitals :math:`\phi_n^\textrm{KS}(\mathbf{r})` with orbital energies :math:`\varepsilon_n^\textrm{KS}$`

.. math::
    :label: equ:KS

    \left\{ -\frac{\hbar^2}{2m}\nabla^2 + V_\mathrm{ext}(\mathbf{r}) + V_\textrm{H}(\mathbf{r}) +V_\textrm{xc}(\mathbf{r})\right\}\phi_n^\textrm{KS}(\mathbf{r})  =\varepsilon_n^\textrm{KS} \phi_n^\textrm{KS}(\mathbf{r}) ,

where :math:`V_\textrm{ext}` is the external potential, :math:`V_\textrm{H}` the Hartree potential, and :math:`V_\textrm{xc}`` the exchange-correlation potential. VOTCA-XTP also contains functionality for projector-based-embedding DFT-in-DFT ground state calculations [Manby:2012]_, in which a chosen _active_ subregion of a molecular system is embedded into an inactive one, reproducing the total energy of the full system ground state exactly.  


Hopping transport
*****************


