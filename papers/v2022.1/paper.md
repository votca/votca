---
title: 'VOTCA'
tags:
  - coarse-graining
  - DFT
  - Quantum calculatons
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


Many physics problems are multi-scale in nature. And the versatile object-oriented
Toolkit for Coarse-grained Applications (`VOTCA`) provides a multi-scale framework for 
liquids and soft-matter system. For instance to understand the exciton dynamics one
needs a whole hierarchy of models reaching from quantum over atomistic and coarse-grained
to reach scales relevant for applications.
For coarse-graining part, VOTCA-CSG, he have added a couple of new methods to derive 
coarse-grained potentials that will be discussed in more details below
For the quantum part, VOTCA-XTP.
Additionally in last two releases we have done major code clean-up and update of the
build systems, but most notable we went to a mono-repository workflow to the user's
and developers' life easier.


# Statement of need

`VOTCA` was originally developed as a platform for method development and comparison.
Since the last software publication in 2015 this stance was strengthened by adding more
methods, more examples and involving more developers.
Many users have used VOTCA to compare different coarse-graining strategies on a neutral
ground and used those insights to use a more specialized package if needed.

(e.g., see
the BOCS10, VOTCA11, DeePCG12, IBIsCO13, MagiC14
,
MSCGFM15 packages).

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

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




