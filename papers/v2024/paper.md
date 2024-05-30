---
title: 'VOTCA: multiscale frameworks for quantum and classical simulations in soft matter'
tags:
  - coarse-graining
  - DFT
  - Quantum calculations
authors:
  - name: Björn Baumeier^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-6077-0467
    affiliation: "1,2"
  - name: Jens Wehner
    orcid: 0000-0003-1203-5000
    affiliation: "1,2,3"
  - name: Nicolas Renaud
    orcid: 0000-0001-9589-2694
    affiliation: 3
  - name: Felipe Zapata Ruiz
    orcid: 0000-0001-8286-677X
    affiliation: 3
  - name: Rene Halver
    orcid: 0000-0002-4895-3762
    affiliation: 4
  - name: Pranav Madhikar
    orcid: 0000-0003-3658-3685
    affiliation: "1,2"
  - name: Ruben Gerritsen
    orcid: 0000-0003-0738-8952
    affiliation: "1,2"
  - name: Gianluca Tirimbo
    orcid: 0000-0002-6641-0761
    affiliation: "1,2"
  - name: Javier Sijen
    affiliation: "1,2"
  - name: David Rosenberger
    orcid: 0000-0001-6620-6499
    affiliation: 5
  - name: Joshua S. Brown
    orcid: 0000-0003-1227-6429
    affiliation: 6
  - name: Vivek Sundaram
    orcid: 0000-0002-9986-8461
    affiliation: "1,2"
  - name: Jakub Krajniak
    orcid: 0000-0001-9372-6975
    affiliation: 7
  - name: Marvin Bernhardt
    orcid: 0000-0001-6520-5208
    affiliation: 8
  - name: Christoph Junghans^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-0925-1458
    affiliation: 9
affiliations:
  - name: Department of Mathematics and Computer Science, Eindhoven University of Technology, the Netherlands
    index: 1
  - name: Institute for Complex Molecular Systems, Eindhoven University of Technology, the Netherlands
    index: 2
  - name: Netherlands eScience Center, the Netherlands
    index: 3
  - name: Forschungszentrum Jülich, Jülich, Germany
    index: 4
  - name: Freie Universität Berlin, Berlin, Germany 
    index: 5
  - name: Oak Ridge National Laboratory, Oak Ridge, TN, USA 
    index: 6
  - name: Independent researcher, Poznań, Poland 
    index: 7
  - name: Technische Universität Darmstadt, Darmstadt, Germany 
    index: 8
  - name: Los Alamos National Laboratory, Los Alamos, New Mexico, USA
    index: 9
  
date: 10 May 2024
bibliography: paper.bib

---

# Summary
Many physical phenomena in liquids and soft matter are multiscale by nature and can involve processes with quantum and classical degrees of freedom occurring over a vast range of length- and timescales. Examples range from structure formation processes of complex polymers or even polymer blends on the classical side to charge and energy transport and conversion processes involving explicit electronic and therefore quantum information. 

The Versatile Object-Oriented Toolkit for Coarse-graining Applications (VOTCA) provides multiscale frameworks built on a comprehensive set of methods for the development of classical coarse-grained potentials (VOTCA-CSG) as well as state-of-the art excited state electronic structure methods based on density-functional and many-body Green's function theories, coupled in mixed quantum-classical models and used in kinetic network models (VOTCA-XTP). 

![Overview of the differnt VOTCA modules and external interfaces. The trajectory reader of VOTCA-CSG in the dashed line box are reused by VOTCA-XTP.\label{fig:VOTCA_graph}](VOTCA_graph.pdf)

# Statement of need

VOTCA was originally developed as a platform for development and comparison of coarse-graining (CSG) methods.
Since the last software publication in 2015 VOTCA-CSG was strengthened by adding more methods, more examples and involving more developers. Many users have used VOTCA to compare different coarse-graining strategies on a neutral ground and, if needed, proceeded with a more specialized package based on the gained insight. Such other coarse-graining packages include, among others, BOCS [@bocs], DeePCG [@deepcg], IBIsCO [@ibisco], MagiC [@magic] and OpenMSCG [@openmscg], some of which are not open-source or specialized in one method. Others have stopped being developed or lack contributions from the greater community. It is difficult to build all-inclusive community package for coarse-graining as it is sometimes hard to consolidate different development styles and constantly changing priorities from sponsors that leave little time for good software engineering practices. In this context we would like to point out that there is a fork of the VOTCA package [@mpip-votca] that contains some feature e.g. Kernel-based machine learning methods [@mlpot], that has not been merged.

Next to strengthening the coarse-graining functionality of VOTCA, another major development direction taken since 2015 is the addition of explicit quantum-mechanical modules aiming at the simulation of static and dynamic properties of electronically excited states in complex molecular environments using multiscale frameworks. Specifically, the VOTCA-XTP part provides an open-source implementation of many-body Green's functions methods (known as $GW$-BSE) with the capability of linking this quantum-mechanical description of a molecular system to a classical (microelectrostatic) environment model and to perform rate-based dynamic simulations with kinetic Monte-Carlo. Basic $GW$-BSE functionality in molecular settings has also more recently been supported in other packages, such as Turbomole [@turbomole], ADF [@ADF2001], PySCF [@pyscf], or MOLGW [@molgw], but these are either commercial or do not provide links to a multiscale framework for complex environments and dynamics.

## Coarse-Graining
In the coarse-graining part of VOTCA, VOTCA-CSG, we made a lot of improvements to the inverse Monte Carlo (IMC) method and have added a new iterative approach the so-called iterative integral equation (IIE) method, which are both described in detail below and in reference therein.

### Inverse Monte Carlo updates
The inverse Monte-Carlo Method introduced by Lyubartsev and Laaksonen in 1995 [@lyubartsev_calculation_1995] is a structure-based coarse graining method, whose goal it is to find an effective pair potential between particles, which reproduces the radial distribution function (RDF) of a reference system (ref) at the coarse grained (CG) resolution. IMC has been part of VOTCA since its first release. In the original implementation the pair potential was determined by iteratively solving a set of linear equations:
\begin{equation}
\label{solve_imc_orig}
(A^TA)\Delta U_{ij} = A^T(g_{ij}^{n-1} - g_{ij}^{ref}),
\end{equation}
where $g_{ij}$ is the RDF between particles $i$ and $j$, $n$ indicates the iteration counter, $\Delta U_{ij}$ is the potential update term, and $A$ and $A^T$ are the Jacobian and its corresponding transpose. The Jacobian $A$ is defined as:
\begin{equation}
\label{eqn_A}
A = \frac{\partial g_{ij}}{\partial U_{ij}},
\end{equation}
where $U_{ij}$ is the pair potential between particles i and j.

Rosenberger _et al._ [@rosenberger_comparison_2016], among others [@toth_iterative_2003; @murtola_coarse-grained_2007; @lyubartsev_systematic_2010], demonstrated that the standard IMC method can suffer from numerical instabilities and/or slow convergence. Therefore, a Tikhonov regularization for IMC has been implemented in VOTCA. This regularization changes the linear equations according to [@rosenberger_comparison_2016]:
\begin{equation}
\label{solve_imc}
(A^TA + \lambda I)\Delta U_{ij} = A^T(g_{ij}^{n-1} - g_{ij}^{ref}),
\end{equation}
where $\lambda$ determines the strength of the regularization and $I$ is the identity matrix.
One can perform a singular value decomposition of the Jacobian $A$ to determine an initial value for $\lambda$[@rosenberger_comparison_2016].
As a rule of thumb $\lambda$ should at least be at the order of the smallest singular values squared.

### Iterative Integral Equation method

The iterative integral equation methods are similar to IMC in that they also aim at reconstructing the RDF of a fine-grained reference system with an effective pair potential.
The main difference is in the construction of the Jacobian, which is approximated in IIE methods from integral equation theory.[@delbaryGeneralizedNewtonIteration2020]
For a molecular fluid, where each molecule is mapped to a single bead, using the Ornstein-Zernicke equation and the hypernetted-chain closure, one arrives at the Jacobian inverse with the form of
\begin{equation}\label{eq:dudg}
A^{-1} = \frac{dU}{dg} = \frac{1}{\beta} \left( 1 - \frac{1}{g} - \mathcal{F}^{-1} \left( \frac{1}{(1 + \rho \hat{h})^2}\right) \mathcal{F} \right) .
\end{equation}
Here, $\hat{h}$ is the Fourier transform of $h = g - 1$ and $\mathcal{F}$ is the Fourier operator.
This approximate Jacobian works well for systems with single-bead molecule representations with convergence as fast as IMC, whereas in the general case, convergence is half as fast as IMC.[@bernhardt_stability_2023]
The costly sampling of the IMC Matrix is not needed, only an RDF which is calculated on twice the range as the potential.[@bernhardt_iterative_2021]

### Constraints

When using the IMC or IIE methods described above to find pair potentials, that best reproduce a reference RDF, one can use the Gauss-Newton algorithm and formulate the problem of finding a potential update $\Delta U$ as a minimization
\begin{equation}\label{eq:GN}
\underset{\Delta U}{\text{arg\,min}} \lvert\lvert \Delta g + A \Delta U \lvert\lvert_2
\end{equation}
where $\Delta g = g - g_\mathrm{target}$.
In that case, additional constraints can be introduced.
For example, it is possible to describe the pressure of a system in terms of the RDF $g$ and the pair potential $U$.
From a target pressure and the current pressure at each iteration, a constraint of the form $B \Delta U = d$ can be described and the constraint is enforced by elimination.[@bernhardt_iterative_2021]



## Electronic Excitations

The most substantial new feature in the VOTCA package is the addition of explicit quantum-mechanical functionalities in the VOTCA-XTP part. The added methods aim at a first-principles-based multiscale modeling of electronically excited states and their dynamics in complex molecular systems. We very briefly describe the three main modules of XTP in the following.

### Density-Functional Theory 
Excited state calculations require a reference ground state calculation within density-functional theory. VOTCA-XTP provides both an automated interface to the ORCA package [@neeseORCAProgramSystem2012] and a lightweight internal DFT engine based on atom-centered Gaussian-type orbitals for method developing and testing. It solves the Kohn-Sham Equations for the molecular orbitals $\phi_n^\textrm{KS}(\mathbf{r})$ with orbital energies $\varepsilon_n^\textrm{KS}$
\begin{equation}
\left\{ -\frac{\hbar^2}{2m}\nabla^2 + V_\mathrm{ext}(\mathbf{r}) + V_\textrm{H}(\mathbf{r}) +V_\textrm{xc}(\mathbf{r})\right\}\phi_n^\textrm{KS}(\mathbf{r})  =\varepsilon_n^\textrm{KS} \phi_n^\textrm{KS}(\mathbf{r}) ,
\label{equ:KS}
\end{equation}
where $V_\textrm{ext}$ is the external potential, $V_\textrm{H}$ the Hartree potential, and $V_\textrm{xc}$ the exchange-correlation potential. VOTCA-XTP also contains functionality for projector-based-embedding DFT-in-DFT ground state calculations [@manby2012simple], in which a chosen _active_ subregion of a molecular system is embedded into an inactive one, reproducing the total energy of the full system ground state exactly.  

### Many-Body Green's Functions and the Bethe-Salpeter Equation
Using the ground-state reference, many-body Green's functions theory with the $GW$ approximation first calculayes _single-particle excitations_ (electron addition or removal) as solutions to the _quasiparticle equations_
\begin{equation}
  \left\{ -\frac{\hbar^2}{2m}\nabla^2 + V_\textrm{ext}(\mathbf{r}) +
    V_\textrm{H}(\mathbf{r})\right\}\phi_n^\textrm{QP}(\mathbf{r}) +
  \int{\Sigma(\mathbf{r},\mathbf{r}',\varepsilon_n^\textrm{QP})\phi_n^\textrm{QP}(\mathbf{r}')d\mathbf{r}'}
  = \varepsilon_n^\textrm{QP} \phi_n^\textrm{QP}(\mathbf{r}) .
\label{equ:QP}
\end{equation}
In place of the exchange-correlation potential in Eq.\ref{equ:KS}, the energy-dependent self-energy operator $\Sigma(\mathbf{r},\mathbf{r}',E)$ occurs in the QP equations. This operator is evaluated using the one-body Green's function in quasi-particle approximation
\begin{equation}
  G(\mathbf{r},\mathbf{r}',\omega) = \sum_n{\frac{\phi_n(\mathbf{r})\phi_n^*(\mathbf{r}')}{\omega-\varepsilon_n+i0^+\textrm{sgn}(\varepsilon_n -\mu)}}
\label{equ:Green}
\end{equation}
as
\begin{equation}
  \Sigma(\mathbf{r},\mathbf{r}',E) = \frac{i}{2\pi} \int{e^{-i\omega 0^+}G(\mathbf{r},\mathbf{r}',E-\omega)W(\mathbf{r},\mathbf{r}',\omega)\,d\omega},
  \label{equ:sigma}
\end{equation}
where $W$ denotes the dynamically screened Coulomb interaction. Assuming that $\phi^\textrm{QP}_n\approx \phi^\textrm{KS}_n$, the quasiparticle energies can be evaluated perturbatively according to
\begin{equation}
  \varepsilon_n^\textrm{QP}= \varepsilon_n^\textrm{KS} + \Delta \varepsilon_n^{GW} =
  \varepsilon_n^\textrm{KS} + \left\langle\phi^\textrm{KS}_n\left\vert
  \Sigma(\varepsilon_n^\textrm{QP})-V_\text{xc} \right\vert\phi^\textrm{KS}_n\right\rangle .
  \label{equ:theory:gw_sc}
\end{equation}
As the correction $\Delta \varepsilon_n^{GW}$ itself depends on $\varepsilon_n^\textrm{QP}$, Eq.\ref{equ:theory:gw_sc} needs to be solved self-consistently.

Neutral excitations with a conserved number of electrons can be obtained from the Bethe-Salpeter Equation (BSE) by expressing coupled electron-hole amplitudes of excitation $S$ in a product basis of single-particle orbitals, i.e., 
\begin{equation}
  \chi_S(\mathbf{r}_\textrm{e},\mathbf{r}_\textrm{h})=\sum_{v}^{\mathrm{occ}}\sum_c^{\mathrm{unocc}}A_{vc}^S\phi_{c}(\mathbf{r}_\textrm{e})\phi^*_{v}(\mathbf{r}_\textrm{h})
  +B_{vc}^S\phi_{v}(\mathbf{r}_\textrm{e})\phi^{*}_{c}(\mathbf{r}_\textrm{h}),
  \label{equ:bsewf}
\end{equation}
where $\mathbf{r}_\textrm{e}$ ($\mathbf{r}_\textrm{h}$) is for the electron (hole) coordinate and $A_{vc}$ ($B_{vc})$ are the expansion coefficients of the excited state wave function in terms of resonant (anti-resonant) transitions between occupied $v$ and unoccupied $c$ states, respectively. In this basis, the BSE turns into an effective two-particle Hamiltonian problem of the form
$$
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
\end{pmatrix}
\label{equ:theory:bseeigenvalue} .
$$

Specifically, the matrix elements of the blocks $\underline{\mathbf{H}}^{\text{res}}$ and $\underline{\mathbf{K}}$ are calculated as
\begin{align}
 H^{\text{res}}_{vc,v'c'}&=D_{vc,v'c'}+\eta K^\mathrm{x}_{vc,v'c'}+K^\mathrm{d}_{vc,v'c'}\label{equ:theory:bseblockH}\\
 K_{cv,v'c'}&=\eta K^\mathrm{x}_{cv,v'c'}+K^\mathrm{d}_{cv,v'c'}\, ,\label{equ:theory:bseblockK}
\end{align}
with 
\begin{align}
D_{vc,v'c'}&=(\varepsilon_c-\varepsilon_v)\delta_{vv'}\delta_{cc'}\label{equ:theory:D},\\
K^\text{x}_{vc,v'c'}&=\iint  \phi_c^*(\mathbf{r}_\textrm{e})\phi_v(\mathbf{r}_\textrm{e})v_{\mathrm{C}}(\mathbf{r}_\textrm{e},\mathbf{r}_\textrm{h}) \phi_{c'}(\mathbf{r}_\textrm{h})\phi_{v'}^*(\mathbf{r}_\textrm{h}) d^3\mathbf{r}_\textrm{e} d^3\mathbf{r}_\textrm{h}\\
K^\text{d}_{vc,v'c'}&=-\iint
               \phi_c^*(\mathbf{r}_\textrm{e})\phi_{c'}(\mathbf{r}_\textrm{e})W(\mathbf{r}_\textrm{e},\mathbf{r}_\textrm{h},\omega=0) \phi_v(\mathbf{r}_\textrm{h})\phi_{v'}^*(\mathbf{r}_\textrm{h})d^3\mathbf{r}_\textrm{e} d^3\mathbf{r}_\textrm{h}
               \, .\label{equ:theory:Kd}
\end{align}
and $\eta=2$ ($\eta=0$) for singlet (triplet) excitations.
Here, $K^\text{x}$ is the repulsive exchange interaction originating from the bare Coulomb term $v_\mathrm{C}$, while the direct interaction $K^\text{d}$ contains the attractive, but screened, interaction $W$ between electron and hole, causing the binding of the electron-hole pair. In Eq.\ref{equ:theory:Kd} it is assumed that the dynamic properties of $W(\omega)$ are negligible, and the computationally less demanding static approximation $\omega=0$ is employed. 


### Quantum-Classical Embedding
Polarization effects of an environment can have significant impact on electronic excitations. As polarization effects are long-ranged accounting for them requires the treatment of large systems which is infeasible with explicit quantum methods such as DFT-$GW$-BSE. Instead, the system is split into a small part with to electronically active subsystem to be treated at quantum (QM) level and a large environment part in which electrostatic and polarization effects are accounted for in classical models (MM). In VOTCA-XTP the QM/MM scheme employs distributed atomic multipole representations for molecules in the MM region, which allows treatment of both the effects of static electric fields and the polarization response as a self-consistent reaction field. Specifically, this classical MM energy for the system is evaluated as 
\begin{equation}
  E_\text{MM} = \frac{1}{2}\sum_{\substack{A,B\\A\neq B}}\sum_{a\in A}\sum_{b\in B}\sum_{tu}{(Q_t^{a} + \Delta Q_t^{a}) T_{tu}^{ab} Q_u^{b}}, 
  \label{equ:EMM}
  \end{equation}
where $A$ and $B$ indicate individual molecules in the system, $a$ and $b$ atoms in the respective molecules, $Q^a_t$ are the static atomic multipole moments of rank $t$ associated to atom $a$, and $T_{tu}^{ab}$ is the tensor describing the interactions between the multipoles moments $Q^{a}_{t}$ and $Q^b_u$ [@stone_distributed_2005]. The induced moments $\Delta Q_t^a$ are generated by the electric field created by moments $t'$ of atom $a' \neq a$ in molecule $A$ and the one generated by the moment $u$ of atom $b$ in molecule $B$:
\begin{equation}
\Delta Q_{t}^{a} = - \sum_{\substack{A,B \in \mathcal{S}\\A\neq B}} \sum_{b \, \in B}\sum_{\substack{a' \in A\\a'\neq a}} \sum_{tt'u}\alpha_{tt'}^{aa'}  T_{t'u}^{a'b} (Q_u^{b} + \Delta Q_u^{b}),
\label{equ:induced}
\end{equation} 
with $\alpha_{tt'}^{aa'}$ the isotropic atomic polarizability on each site. To avoid effects of spurious overpolarization, a damped version of the interaction tensor (Thole damping [@stone_distributed_2005]) is used. Then, the static and induced multipoles in the MM region also interact with the electron density in QM region via an additional external potential to Eq.\ref{equ:KS}. At the same time, the explicit electrostatic field from the QM density is included in polarizing the MM region.

# Code Structure

For the last couple of years, we have also focused on code hardening and the introduction of better software engineering practices. Original VOTCA was designed as modules in separate repositories, but as many other projects, this turned out to be quite cumbersome hence we switched to a mono-repo. With recent performance improvements in the git tools, the benefits of a single repository by far out-weigh the downside of the very complex workflow of multiple repositories. The module structure still exists in the source code.

Additionally, we have added continuous integration testing through GitHub action for 50+ different compiler and operating system combinations. The also perform continous deployment to the GitHub Docker registry.
And releases get rolled into all major linux distributions, HomeBrew, Spack and FreeBSD.

## Code Modernization

We did a lot of code refactoring and bumped the C++ standard to 17. We also modernized our usage of CMake and switched to a mostly target-base scheme. An attempt to port our particle structure on top of Cabana [@cabana] was made, due to incompatibilties between Kokkos [@kokkos] and Eigen, we will delay this effort.

## Updates in VOTCA-CSG

The particle and molecule data structure were refactored, and we add support of the H5MD format, which is described below in details.
 
### H5MD support

The recent version of VOTCA supports the `H5MD` [@debuyl2014h5md] file format, which internally uses `HDF5` [@hdf5] storage. This is a very fast and scalable method for storing molecular trajectories, already implemented in simulation packages such as `LAMMPS`, `ESPResSo++`, and `ESPResSo`. VOTCA recognizes the trajectory file format by the extension. In the case of H5MD, it expects a `.h5` extension. Following the H5MD concepts, the particle trajectories are organized in the `particles` container. This container can handle multiple subsets of the studied system. Therefore, we must define `h5md_particle_group` in the XML topology file to declare which subset of particles to use. The reader handles both coordinates (if present), forces, and velocities.

## Design of VOTCA-XTP
Data structures related to atomistic properties (topology, molecules, segments, fragments, atoms) in XTP are reused or build upon those of CSG. Linear algebra related structures and functionalities are handled by Eigen [@eigenweb] which can be accelerated by internally calling the Intel Math Kernel Library [@intelMKL]. Exchange-correlation functionals are provided by the Library of eXchange-Correlation (`LIBXC`) functionals [@Lehtola2018], while `libint` [@Libint2] and `libecpint` [@Shaw2017] are used for the evaluation of molecular integrals of many-body operators over Gaussian functions.

VOTCA–XTP is designed as a library, which is linked to very thin executables. These executables provide a variety of calculators by adding keywords on the command line. Virtual interfaces and factory patterns make the addition of new calculators simple. The same architecture is used for external DFT and MD codes, making VOTCA–XTP easily extensible. Lower-level data structures make use of template metaprogramming to support a variety of data types. VOTCA-XTP provides different functionalities in three types of _calculator_ classes: 

- a collection of tools that do not require information of a mapped MD trajectory, including a specific DFT-$GW$-BSE calculator in `tools` callable by `xtp_tools` 
- analysis and not-high-throughput applications that require a mapped MD trajectory in `calculators` callable by `xtp_run`
- high-throughput, high-performance applications that require a mapped MD trajectory in `jobcalculators` callable by `xtp_parallel`

In general, VOTCA-XTP uses shared-memory parallelization in the heavy calculations involving the quantum methods, with the possibility to seamlessly offload matrix-matrix and matrix-vector operations to GPU via `CUDA`. 
  
### PyXTP

The PyXTP python package distributed with VOTCA, contains python bindings to the main functionalities of VOTCA-XTP. These python bindings were created using pybind11 [@pybind11] and provide users with a dedicated Atomistic Simulation Environment (ASE) [@ase] calculator. The use of this calculator not only facilitates the adoption of VOTCA-XTP by non-experts users, but they also allow integrating VOTCA-XTP in the broader ASE ecosystem. 

The following snippet of code illustrate the use of PyXTP. This small code optimize the geometry of a CO molecule in the first excited singlet states. As seen in the code, the XTP calculator is used to compute the forces on the nuclei while the geometry optimization itself is driven by ASE functionalities. 

```python
from pyxtp import xtp
from ase.io import write
from ase.build import molecule
from ase.optimize import QuasiNewton

# create a distorted CO molecule
atoms = molecule('CO')
atoms.rattle()

# instantiate the calculator
calc = xtp(nthreads=2)

# select the state for which to compute the forces 
calc.select_force(energy='singlets', level=0, dynamic=False)

# this allows to change all options
calc.options.dftpackage.functional = 'PBE'
calc.options.dftpackage.basisset = 'def2-svp'
calc.options.dftpackage.auxbasisset = 'aux-def2-svp'

# set up the logger
calc.options.logging_file = 'CO_forces.log'

# set the calculator
atoms.calc = calc

# optimize the geometry
dyn = QuasiNewton(atoms, trajectory='test.traj')
dyn.run(fmax=0.01)
write('final.xyz', atoms)
```

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project. 
We acknowledge support by the Innovational Research Incentives Scheme Vidi of the Netherlands Organisation for Scientific Research (NWO) with project number 723.016.002. Funding is also provided by NWO and the Netherlands eScience Center for funding through project number 027.017.G15, within the Joint CSER and eScience program for Energy Research (JCER 2017).
Los Alamos National Laboratory (LANL) is operated by Triad National Security, LLC, for the National Nuclear Security Administration of US Department of Energy (Contract No. 89233218CNA000001).
This work is approved for unlimited release with report number LA-UR-24-25313.

# References




