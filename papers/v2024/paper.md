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

## Electronic Excitations
**Bjoern**


# Code Structure
**Some general statements**
## CSG
## XTP
- builds on CSG
- calculator concept
- tools
- parallelism
- GPU (**ask eScience? Jens?**)
## PyXTP & ASE
**ask Nico to write a bit**
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




