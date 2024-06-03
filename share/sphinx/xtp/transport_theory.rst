Charge and Energy Transport Model
#################################

Partitioning of the system
**************************

Before the charge or exciton dynamics can be calculated using rate models, a given morphology has to be divided into sites, on which excitations localize. This requires knowledge of the localization of the electronic states, as each site should preferably contain a single diabatic state of the complete system. In detail knowledge of the localization requires electronic structure calculations of the whole system, which are not feasible for larger systems. Instead, the partitioning into segments can only be justified a posteriori by assessing the magnitude of the electronic couplings between segments. 

.. _fig_theory_fragment_def:
.. figure:: fig/fragment_def.svg
    :width: 400
    :align: center

    Partitioning scheme as employed in VOTCA. a) For small, inflexible molecules, the whole molecule is replaced by a rigid fragment. b) For larger molecules, more fragments are used to allow for some flexibility for the slower degrees of freedom. The whole molecule is still one segment. c) For a polymer, the chain is cut into segments, using input from quantum mechanical studies. Each segment may also contain more than one rigid fragment.

In many systems chemical intuition allows defining partitioning, which satisfies this test retrospectively, e.g., in most organic electronic systems a conjugated system is a good minimal guess. As shown in :numref:`fig_theory_fragment_def` all atoms belonging to a single site are collectively referred to as a conjugate segment. 

In small molecule systems each molecule is typically chosen as a single conjugated segment, because electronic structure calculations typically show complete delocalization of electronic states over the individual molecules. In polymeric systems this simple approach may fail, as torsion of the backbone and electron-phonon interaction may lead to localization on polymeric subunits instead of delocalization over the whole polymer chain. In these cases quantum mechanical calculations of the whole chain are necessary to determine localization sites.

After the system is partitioned into segments, these segments are additionally separated into rigid fragments, as depicted in :numref:`fig_theory_fragment_def`. These fragments are typically insensitive to thermal fluctuations, e.g., a thiophene or benzene ring stays planar even at elevated temperatures. Although the fragment definitions are theoretically not necessary for the calculation of rates, the molecular dynamics configurations of the rigid fragments are often replaced with quantum mechanically optimized rigid fragments, in order to *integrate out* fast degrees of freedom, e.g. the promoting modes, which are already described in the rate models. This has the additional benefit of removing some of the mismatch between quantum mechanically obtained geometries and molecular dynamics results.  


Site Energies
*************

The site energy, :math:`\Delta E_a`, of a molecule :math:`a` is basically defined via :eq:`equ:coarsesiteE`. Although this is a simple definition it is frayed with computational difficulties as the energy :math:`U_a` not only depends on the isolated molecule but also the chemical environment it is located in. As discussed before, the electronic state also must be sufficiently localized on the molecule to even make the definition of a site energy valid. Sufficient localization entails weaker coupling to the surrounding, allowing us to treat the chemical environment as a small perturbation to the electronic state. 

Although there are exchange, dispersion, and electrostatic inrteractions between molecules, we focues only on the last, electrostatic interaction. There are a couple of reasons for this approximation. Firstly, we are not really interested in the absolute values of site energies but in differences between different electronic states, e.g. :math:`\Delta E_a=U^1_a-U^0_a`. Consequently, all interactions, which influence both states nearly equally, e.g. dispersion interaction, cancel out. Secondly, electrostatic interaction is very directional, thus it is very sensitive to changes in the charge distribution on the molecules. Thirdly, and most importantly it is typically the strongest interaction present, due to its slow decay with distance. So the site energy of a molecule :math:`a` in state :math:`i` is decomposed into the energy from the internal degrees of freedom :math:`E_\text{int}` and the electrostatic contribution from the environment :math:`E_\text{estat}`:

.. math::
    :label: equ:coarse:siteenergydiff

    U_a^i(\xi_a^i)\approx E_\text{int}^i+E_\text{estat}^i

The internal site energy difference :math:`\Delta E_{\text{int}}=E_{\text{int}}^1-E_{\text{int}}^0` can be calculated from the vacuum PES of the single molecule via:

.. math::
    :label: equ:interalsiteenergy

    \Delta E_{\text{int}}= U^1_\text{vac}(\xi^1_a)-U^0_\text{vac}(\xi^0_a).

Theoretically :math:`\Delta E_a` is an adiabatic excitation energy calculated from the respective geometries for the ground and excited state. Depending on the electronic excitation (electron, hole, exciton) and the system at hand, slow (compared to the timescale of hopping) vibrational modes may seriously change the local excitation energy by distorting molecular geometries. These distortions lead to disorder in the internal site energies. Using vacuum optimized geometries may not be a good approximation, as they optimize out this dynamic, but on the timescale of electronic motion, frozen disorder. In practice geometries from molecular dynamics simulations are used for the ground and excited states, using the mapping procedure described earlier to account for the geometry reorganization upon excitation in the site energies.

As mentioned above, the electrostatic contribution to the site energy of a molecule :math:`E_{\text{estat}}^i` depends on the surrounding. The question is how to treat this effect. A fully quantum mechanical treatment of the whole system is impossible. Instead, the interaction between the quantum mechanical excitation and the environment can be modeled on different levels (see :numref:`fig_theory_siteenergies`). Starting with the quantum mechanical excitation in vacuum, :numref:`fig_theory_siteenergies` (a), either the electrostatic potential of the surrounding molecules can be fully added to the Hamiltonian operator of the quantum mechanical region :numref:`fig_theory_siteenergies` (c). Instead the interaction can be described as classically interacting charge densities, :numref:`fig_theory_siteenergies` (b), which will be described in the next paragraphs. The effects of polarization contributions, :numref:`fig_theory_siteenergies` (d), which allow the environment to react to excitations will be described later.

.. _fig_theory_siteenergies:
.. figure:: fig/siteenergies.svg
    :width: 400
    :align: center

    Overview of the site energy calculation options in the VOTCA suite using either a perturbative classical approach or a coupled QM/MM approach. For the classical approach the initial starting point is a vacuum QM calculation, yielding the ground state and excited state energies in vacuum (a). Adding purely classical electrostatic leads to different environments for every molecule and thus a spread of the HOMO and LUMO energies (b). (d) Adding classical polarization (indicated by black arrows) on neighboring molecules allows the system to relax and reduce the energetic disorder. Likewise, the QM/MM treatment inserts a fully quantum mechanical description of the molecule surrounded by classical charges (c) or as in the classical model a full polarization response of the environment (d). The difference being, that a the polarization response of the molecule in question is more flexible in the QM/MM treatment, as charge can flow over the whole molecule.

Focusing on one molecule, the electrostatic contribution of the environment to its Hamiltonian reads:

.. math::
    :label: equ:theory:extham
 
    \hat{H}_\text{ext}=\int d^3{\mathbf{r}} \hat{\rho}(\mathbf{r})\varphi(\mathbf{r})
 
where :math:`\hat{\rho}(\mathbf{r})` is the molecules charge density operator and :math:`\varphi(\mathbf{r})` is the external electrostatic potential due to the environment. It can be written as: 

.. math::
    :label: equ:theory:extpot
 
    \varphi(\mathbf{r})=\sum_K \int d^3{\mathbf{r'}}\frac{\rho_K(\mathbf{r'})}{|\mathbf{r'}-\mathbf{r}|},

where :math:`K` runs over all molecules in the environment and :math:`\rho_K(\mathbf{r'})` is the molecular charge density belonging to molecule :math:`K`. Here, we modeled the neighboring molecules as classical entities by replacing their charge density operator with the respective classical charge densities.

Calculating the external potential, :eq:`equ:theory:extpot`, requires an integration over the charge densities of the adjacent molecules, which although possible is computationally very demanding for large systems. Instead, a distributed multipole expansion is used to calculate the external potential. In a distributed multipole expansion, the molecular charge density :math:`\rho(\mathbf{r})` is not expanded via multipoles centered around a single point. Instead, the molecule is divided into regions, which are described by their own multipole moments. A region typically comprises a single atom, with the origin located at the atom. The potential of a single region of molecule can be expanded as:

.. math::

    \varphi(\mathbf{r})^\text{el}=T(\mathbf{r})q-T_\alpha(\mathbf{r})d_\alpha +\frac{1}{3}T_{\alpha\beta}(\mathbf{r})\Theta_{\alpha\beta}-\dots+\frac{(-1)^n}{(2n-1)!!}T_{\alpha\beta\dots\nu}^{(n)}(\mathbf{r})\zeta^{(n)}_{\alpha\beta\dots\nu},

with :math:`q` being the charge of that region, :math:`d` the dipole moment, :math:`\Theta_{\alpha\beta}` the Cartesian quadrupole tensor and :math:`\zeta^{(n)}_{\alpha\beta\cdots\nu}` the corresponding higher moments. The Cartesian interaction tensor :math:`T_{\alpha\beta\cdots\nu}^{(n)}` is defined as:

.. math::

    T_{\alpha\beta\cdots\nu}^{(n)}(\mathbf{r})=\nabla_\alpha\nabla_\beta\dots\nabla_\nu\frac{1}{|\mathbf{r}|}



Finally, :math:`\varphi(\mathbf{r})` reduces to a sum over all regions :math:`r(K)` of all adjacent molecules :math:`K`:

.. math::

    \varphi(\mathbf{r})=\sum_K\sum_{r(K)} \varphi(\mathbf{r})^\text{single}_{K,r(K)}


The resulting potential can then be inserted into :eq:`equ:theory:extham` and in conjunction with the *ab-initio* methods the excited and ground state energy can be calculated. This corresponds to the first QM/MM approach in :numref:`fig_theory_siteenergies`, where the quantum mechanical charge density reacts to the static environment. This approach still has two drawbacks. Firstly, performing an *ab-initio* calculation for every molecule in a large system can be computationally infeasible. Secondly, so far we have assumed that the environment does not react to the excitation. In reality, the environment tries to lower the total energy by screening the excitation, e.g. the environment polarizes. The polarization response of the environment in turn will lead to a polarization of the center molecule. This feedback loop requires a self-consistent treatment, regardless of whether the central molecule is described classically or quantum mechanically.

To address the first problem, we simply substitute the quantum mechanical treatment of the center molecule by a classical description via distributed multipole moments. The electrostatic interaction :math:`V^{AB}_\text{elstatic}` between two molecules :math:`A` and :math:`B` then reduces to:

.. math::
    :label: equ:theory:classicinteraction

    \begin{align}
    V^{AB}_\text{elstatic}=&\sum_{a\in A}\sum_{b\in B}\left[T^{ab}q^aq^b+T^{ab}_\alpha(q^ad^b_\alpha-d_\alpha^aq^b)+\right.\nonumber\\
    &\left.T^{ab}_{\alpha\beta}(\frac{1}{3}q^a\Theta^b_{\alpha\beta}-d^a_{\alpha}d^b_{\beta}+\frac{1}{3}\Theta^a_{\alpha\beta}q^b)+\dots\right]
    \end{align}

with :math:`T^{ab}=T(\mathbf{r}_a-\mathbf{r}_b)`.

The electrostatic contribution to the site energy of molecule :math:`A` in state :math:`i` then reads:

.. math::
    :label: equ:coarse:electrostatic

    E^i_{\text{estat}}=\sum_{B}V^{AB}_\text{elstatic}(a^i),

where :math:`(a^i)` denotes that the multipoles belonging to the state :math:`i` should be used. :math:`B` runs over all other molecules, with their respective ground state multipoles. The classical environment correction for a transition from the ground state to :math:`n`-th excited state then be :math:`\Delta E=E^0_{\text{estat}}-E^n_{\text{estat}}`.

The second problem of how to incorporate the polarization response of the environment requires an extension of the distributed multipole model.
Every atom :math:`a` of a molecule is also assigned a polarization tensor :math:`\alpha^a_{\alpha\beta}`. An external field :math:`F^a_\alpha` will lead to an induced dipole moment on site :math:`a`. The induced dipoles are not only caused by the external field originating from the by static multipoles on other atoms but also by the fields originating from induced dipole moments on other atoms. The total field at :math:`a` thus is a combination of the static external field :math:`F^a_\alpha` as well as the induced fields from the induced dipoles on all other atoms. The induced dipole on :math:`a` is: 

.. math::
    :label: equ:theory:apple

    \Delta d^a_\gamma=\alpha^a_{\alpha\beta}\left(F^a_\alpha+\sum_{b\neq a} T^{ab}_{\alpha\beta}\Delta d^b_\beta \right)

The assignment of atomic polarizabilities is called the Applequist model~\cite{Applequist_Atomdipoleinteraction_1972}. Here, the polarizations are only correlated via the fields they generate (see :eq:`equ:theory:apple`) and cannot model the flow of charges across the whole molecule. For molecules with large conjugated systems, which can displace charges over larger distances, necessitating more advanced polarization models or a full quantum mechanical treatment are necessary. The system of equations, :eq:`equ:theory:apple`, has to be solved self-consistently or via matrix inversion, which both increase the computational cost approximately by one order of magnitude in comparison to a static approach.

Letting the environment react to the excitation on the center molecule explicitly couples the QM and MM system to each other. So energy differences cannot be calculated for the embedded molecule but instead total energies have to be subtracted. For example, the first excitation energy :math:`\Omega_1` has to be calculated via:

.. math::
    :label: equ:theory:qmmmexcitation

     \Omega_1=\left[GS(BG[\Omega_1])+\Omega_1(BG[\Omega_1])\right]-GS(BG[GS])

with :math:`GS(BG[\Omega_1])` being the ground state energy, in a background/environment polarized for the first excited state. As a result, for each site energy the excited state in its polarization environment as well as the ground state in its polarization environment have to be calculated, further increasing the computational cost in comparison to a simple static approach.

Finally, the induced-induced interaction of the Applequist model leads to a polarization catastrophe at short distances. This is an artifact of the point-dipole description of polarization, as in reality the polarization response is smeared out. An empirical fix, called Thole damping, avoids this overpolarization by introducing a smeared out charge density :math:`\varrho(u)` with :math:`u=r/S_{ab}` and :math:`S_{ab}=(\alpha_a\alpha_b)^{1/6}`, which modifies the dipole-dipole interaction between atoms :math:`a` and :math:`b` with mean polarizabilities :math:`\alpha_a`, :math:`\alpha_b` as follows:

.. math::

    \tilde{T}^{ab}_{\alpha\beta}=\frac{f_e R^2\delta_{\alpha\beta}-3f_tR_\alpha R_\beta}{R^5}.

Here :math:`f_e` and :math:`f_t` are the corresponding damping functions, which depend on the form of the smeared out charge density. We use the smeared density from the AMOEBA force field

.. math::

    \varrho(u)=\frac{3a}{4\pi}\exp(-au^3).

This leads to :math:`f_e=1-\exp(-\lambda u^3)` and :math:`f_t=1-(1+\lambda u^3)\exp(-\lambda u^3)`. In this work we use :math:`\lambda=0.37`. The atomic polarizabilities :math:`\alpha` were taken from a generic set of element-based polarizabilities. These are sufficient to reproduce the molecular polarizability of many organic compounds. However, :math:`\pi`-conjugated molecules with their delocalized :math:`\pi`-orbitals have higher polarizabilities which cannot be reproduced with the generic atomic polarizabilities. In this case the atomic polarizabilities are scaled to match the molecular polarizabilities from ab-initio calculations for the respective compounds.  

A technical difficulty that may even arise from the simple electrostatic summation in :eq:`equ:theory:classicinteraction` is the question of how many neighboring molecules have to be considered. Traditionally a cut-off based approach is used, taking all molecules within a certain distance into account. This is problematic because the electrostatic interaction does not converge with increasing cutoff size. Although this is ameliorated by the fact, that most systems are not charged and thus the charge-charge interaction cancels out, even charge-quadrupole interaction is only conditionally convergent. Depending on the dimensionality of the system and the excitation of interest more advanced methods like *long-range aperiodic embedding* have to be used.

This finally leaves the question of how to reduce the molecular density :math:`\rho(r)` to a set of multipole moments located at the atoms :math:`\{Q_{\alpha\beta}^i\}`. There is a number of different methods to arrive at atom-centered multipoles, from Mulliken charges [Mulliken:1955]_ over [Reed:1985]_ charges to GDMA [Stone:2005]_ multipoles. We will briefly outline the CHELPG (CHarges from ELectrostatic Potential, Grid-base) method [Chirlian:1987]_ , which, although it does not produce chemically intuitive charges, is designed to reproduce the electrostatic potential outside the molecule optimally via point charges. The point charges derived from molecular charge densities are typically referred to as partial charges.

CHELPG partial charges :math:`\{q_i\}` are derived by calculating the electrostatic potential :math:`\phi_{el}(\mathbf{r})` from the electron density :math:`\rho(\mathbf{r})` at a large number of grid points outside the molecule :math:`\{\mathbf{g}_i\}` and then adjusting the partial charges to reproduce the potential in a least square sense as optimally as possible, e.g. minimizing:

.. math::

    \text{min}=\left(\sum_{j=0}^K \phi_{el}(\mathbf{g}_j)-\sum_{i=0}^N \frac{q_i}{|\mathbf{r}_i-\mathbf{g}_j|}\right)^2-\lambda \left(q_{\text{mol}}-\sum_{i=0}^N q_i\right),


where :math:`j` runs over all grid points, :math:`\phi_{el}(\mathbf{g}_j)` is the respective potential at that grid point and :math:`N` is the number of atomic sites. :math:`\lambda` is a Lagrange multiplier to constrain the optimization to the desired total charge of the molecule :math:`q_{\text{mol}}`.