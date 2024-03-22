---
title: 'VOTCA: multiscale frameworks for quantum and classical simulations in soft matter'
tags:
  - coarse-graining
  - DFT
  - Quantum calculations
authors:
  - name: Bjoern Baumeier # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-6077-0467
    affiliation: "1, 2"
  - name: Christoph Junghans^[corresponding author] # note this makes a footnote saying 'co-first author'
    affiliation: 3
  - name: Nicolas Renaud
    affiliation: 4
  - name: Felipe Zapata Ruiz
    affiliation: 4
  - name: Rene Halver
  - name: Pranav Madhikar
    affiliation: "1,2"
  - name: Ruben Gerritsen
    affiliation: "1,2"
  - name: Gianluca Tirimbo
    affiliation: "1, 2"
  - name: Javier Sijen
    affiliation: "1, 2"
  - name: David Rosenberger
  - name: Inti Pelupessy
    affiliation: 4
  - name: Joshua S. Brown
    affiliation: 2
  - name: Vivek Sundaram
  - name: Jakub Kraijniak
  - name: Marvin Bernhardt
affiliations:
  - name: Department of Mathematics and Computer Science, Eindhoven University of Technology, the Netherlands
    index: 1
  - name: Institute for Complex Molecular Systems, Eindhoven University of Technology, the Netherlands
    index: 2
  - name: Los Alamos National Laboratory
    index: 3
date: 23 March 2023
bibliography: paper.bib

---

# Summary
Many physical phenomena in liquids and soft matter are multiscale by nature and can involve processes with quantum and classical degrees of freedom occurring over a vast range of length- and timescales. Examples range from structure formation processes of complex polymers or even polymer blends **(better example linked to classical coarse-graining?)** on the classical side to charge and energy transport and conversion processes involving explicit electronic and therefore quantum information. 

The Versatile Object-Oriented Toolkit for Coarse-grained Applications (`VOTCA`) provides multiscale frameworks built on a comprehensive set of methods for the development of classical coarse-grained potentials (`VOTCA-CSG`) as well as state-of-the art excited state electronic structure methods based on density-functional and many-body Green's function theories, coupled in mixed quantum-classical models and used in kinetic network models (`VOTCA-XTP`). 

### summarize some (technical) details here
- initial release...
- **C+P from JCP, probably too detailed and can go to Code Structure?** VOTCA–XTP is written in C++ and mostly adheres to the C++14 standard.78 It can be obtained on www.github.com/votca/xtp. We use the Git feature branch workflow combined with code review and continuous integration, which executes code formatting, static analyzers, debug and release builds, and the test-suite. We use CMake as the build system, which also resolves the inclusion of external dependencies. The linear algebra is handled by Eigen,79 which can be accelerated by internally calling the Intel Math Kernel Library.80 For serialization, the HDF5 format is used via the canonical libraries.81 Exchange-correlation functionals are provided by the Library of eXchange-Correlation (LIBXC) functionals.82 Various boost packages83 are used for file system string operations. Doxygen is used to document the APIs of VOTCA–XTP and automatically deploys to http://doc.votca.org.
VOTCA–XTP is designed as a library, which is linked to very thin executables, which can execute a variety of calculators by adding keywords on the command line. Virtual interfaces and factory patterns make the addition of new calculators simple. The same architecture is used for external DFT and MD codes, making VOTCA–XTP easily extensible. Lower-level data structures make use of template metaprogramming to support a variety of data types.
- summary of recent CSG developments (added methods in CSG)
- summary of features for XTP part (all new? GTO implementation of GW-BSE, QMMM modes, kinetic MC with finite lifetime, subsystem-embedding)
- Additionally in last two releases we have done major code clean-up and update of the build systems, but most notable we went to a mono-repository workflow to the user's and developers' life easier.


# Statement of need (incl. Capabilities)

`VOTCA` was originally developed as a platform for method development and comparison.
Since the last software publication in 2015 this stance was strengthened by adding more
methods, more examples and involving more developers.
Many users have used VOTCA to compare different coarse-graining strategies on a neutral
ground and used those insights to use a more specialized package if needed.

(e.g., see
the BOCS10, VOTCA11, DeePCG12, IBIsCO13, MagiC14
,
MSCGFM15 packages).

## Coarse-Graining
**Christoph**

### IMC updates - David
The inverse Monte-Carlo Method (IMC) introduced by Lyubartsev and Laaksonen in 1995\cite{lyubartsev_calculation_1995} is a structure-based coarse graining method, whose goal it is to find  an effective pair potential between particles, which reproduces the radial distribution function (RDF) of a reference system (ref) at the coarse grained (CG) resolution. IMC has been part of VOTCA since its first release. In the original implementation the pair potential was determined by iteratively solving a set of linear equations:
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

Rosenberger et al.\cite{rosenberger_comparison_2016}, among others\cite{toth_iterative_2003,murtola_coarse-grained_2007,lyubartsev_systematic_2010}, demonstrated that the standard IMC method can suffer from numerical instabilities and/or slow convergence. Therefore, a Tikhonov regularization for IMC has been implemented in VOTCA. This regularization changes the linear equations according to\cite{rosenberger_comparison_2016}:
\begin{equation}
\label{solve_imc}
(A^TA + \lambda I)\Delta U_{ij} = A^T(g_{ij}^{n-1} - g_{ij}^{ref}),
\end{equation}
where $\lambda$ determines the strength of the regularization and $I$ is the identity matrix.
One can perform a singular value decomposition of the Jacobian $A$ to determine an initial value for $\lambda$\cite{rosenberger_comparison_2016}.
As a rule of thumb $\lambda$ should at least be at the order of the smallest singular values squared.

### IIE method

The iterative integral equation (IIE) methods are similar to IMC in that they also aim at reconstructing the RDF of a fine-grained reference system with an effective pair potential.
The main difference is in the construction of the Jacobian, which is approximated in IIE methods from integral equation theory.\cite{delbaryGeneralizedNewtonIteration2020}
For a molecular fluid, where each molecule is mapped to a single bead, using the Ornstein-Zernicke equation and the hypernetted-chain closure, one arrives at the Jacobian inverse with the form of
\begin{equation}\label{eq:dudg}
A^{-1} = \frac{dU}{dg} = \frac{1}{\beta} \left( 1 - \frac{1}{g} - \mathcal{F}^{-1} \left( \frac{1}{(1 + \rho \hat{h})^2}\right) \mathcal{F} \right) .
\end{equation}
Here, $\hat{h}$ is the Fourier transform of $h = g - 1$ and $\mathcal{F}$ is the Fourier operator.
This approximate Jacobian works well for systems with single-bead molecule representations with convergence as fast as IMC, whereas in the general case, convergence is half as fast as IMC.\cite{bernhardt_stability_2023}
The costly sampling of the IMC Matrix is not needed, only an RDF which is calculated on twice the range as the potential.\cite{bernhardt_iterative_2021}

### Constraints

When using the IMC or IIE methods described above to find pair potentials, that best reproduce a reference RDF, one can use the Gauss-Newton algorithm and formulate the problem of finding a potential update $\Delta U$ as a minimization
\begin{equation}\label{eq:GN}
\underset{\Delta U}{\text{arg\,min}} \lvert\lvert \Delta g + A \Delta U \lvert\lvert_2
\end{equation}
where $\Delta g = g - g_\mathrm{target}$.
In that case, additional constraints can be introduced.
For example, it is possible to describe the pressure of a system in terms of the RDF $g$ and the pair potential $U$.
From a target pressure and the current pressure at each iteration, a constraint of the form $B \Delta U = d$ can be described and the constraint is enforced by elimination.\cite{bernhardt_iterative_2021}







## Electronic Excitations
**Bjoern**


# Code Structure
**Some general statements**
## CSG

### Code Refactor - Josh
### H5MD support

The recent version of VOTCA supports the H5MD\cite{debuyl2014h5md} file format, which internally uses HDF5\cite{hdf5} storage. This is a very fast and scalable method for storing molecular trajectories, already implemented in simulation packages such as LAMMPS, ESPResSo++, and ESPResSo.
VOTCA recognizes the trajectory file format by the extension. In the case of H5MD, it expects a `.h5` extension. Following the H5MD concepts, the particle trajectories are organized in the `particles` container.
This container can handle multiple subsets of the studied system. Therefore, we must define `h5md_particle_group` in the XML topology file to declare which subset of particles to use.
The reader handles both coordinates (if present), forces, and velocities.

## XTP
- builds on CSG
- calculator concept
- tools
- parallelism
- GPU (**ask eScience? Jens?**)
  
### PyXTP

The PyXTP pyhon package distributed with VOTCA, contains python bindings to the main functionalities of `VOTCA-XTP`. These python bindings were created using pybind11 [@pybind11] and provide users with a dedicated Atomisic Simulation Environment (ASE) [@ase] calculator. The use of this calculator not only facilitates the adoption of `VOTCA-XTP` by non-experts users but they also allow to integrate `VOTCA-XTP` in the broader ASE ecosystem. 

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

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures


# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References




