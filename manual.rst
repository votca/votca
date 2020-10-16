Disclamer
=========

This manual is not complete. The best way to start using the software is
to look at provided tutorials. The reference section is generated
automatically from the source code, so please make sure that your
software and manual versions match.

Citations
=========

| Development of this software depends on academic research grants. If
  you are using the package, please cite the following papers
| :raw-latex:`\cite{mashayakrelative}` Relative entropy and
  optimization-driven coarse-graining methods in VOTCA,
| S.Y. Mashayak, Mara Jochum, Konstantin Koschke, N.R. Aluru, Victor
  Rühle, and Christoph Junghans,

| :raw-latex:`\cite{ruhle2011hybrid}` Hybrid approaches to
  coarse-graining using the VOTCA package: liquid hexane,
| Victor Rühle and Christoph Junghans,

| :raw-latex:`\cite{Ruehle:2009.a}` Versatile Object-oriented Toolkit
  for Coarse-graining Applications
| Victor Rühle, Christoph Junghans, Alexander Lukyanov, Kurt Kremer, and
  Denis Andrienko

Development
===========

The core development is currently taking place at the Los Alamos
National Laboratory and Max Planck Institute for Polymer Research,
Mainz, Germany.

Copyright
=========

 votcais free software. The entire package is available under the Apache
License. For details, check the LICENSE file in the source code. The
votcasource code is available on our homepage, .

Introduction
============

Versatile Object-oriented Toolkit for Coarse-graining Applications, or
votca, is a package which helps to systematically coarse-grain various
systems :raw-latex:`\cite{Ruehle:2009.a}`. This includes deriving the
coarse-grained potentials, assessing their quality, preparing input
files required for coarse-grained simulations, and analyzing the latter.

A typical coarse-graining workflow includes *sampling* of the system of
interest, *analysis* of the trajectory using a specific *mapping* and a
coarse-graining *method* to derive coarse-grained potentials and, in
case of iterative methods, running coarse-grained simulations and
iteratively *refining* the coarse-grained potentials.

In most cases, coarse-graining requires canonical sampling of a
reference (high resolution) system. In addition, iterative methods
require canonical sampling of the coarse-grained system. The sampling
can be done using either molecular dynamics (MD), stochastic dynamics
(SD), or Monte Carlo (MC) techniques. The latter are implemented in many
standard simulation packages. Rather than implementing its own MD/SD/MC
modules, votcaallows swift and flexible integration of existing programs
in such a way that sampling is performed by the program of choice. At
the moment, an interface to GROMACS :raw-latex:`\cite{gromacs4}`
simulation package is provided. The rest of the analysis needed for
systematic coarse-graining is done using the package tools.

ht4cm |image|

The workflow can be exemplified on coarse-graining of a propane liquid.
A single molecule of propane contains three carbon and eight hydrogen
atoms. A united atom coarse-grained representation of a propane molecule
has three beads and two bead types, A and B, with three and two
hydrogens combined with the corresponding atom, as shown in
fig. [fig:intro:propane]. This representation defines the , as well as
the bonded coarse-grained degrees of freedom, such as the bond :math:`b`
and the bond angle :math:`\theta`. Apart from the bonded interactions,
:math:`u_b` and :math:`u_\theta`, beads belonging to different molecules
have non-bonded interactions, :math:`u_\text{AA}`, :math:`u_\text{AB}`,
:math:`u_\text{BB}`. The task of coarse-graining is then to derive a
potential energy surface :math:`u` which is a function of all
coarse-grained degrees of freedom. Note that, while the atomistic bond
and angle potentials are often chosen to be simple harmonic functions,
the coarse-grained potentials cannot be expressed in terms of simple
analytic functions. Instead, tabulated functions are normally used.

The coarse-graining *method* defines criteria according to which the
potential energy surface is constructed. For example, for the bond
:math:`b` and the angle :math:`\theta`  can be used. In this case a
coarse-grained potential will be a potential of mean force. For the
non-bonded degrees of freedom, the package provides or methods. In this
case the radial distribution functions of the coarse-grained model will
match those of the atomistic model. Alternatively, (or multiscale
coarse-graining) can be used, in which case the coarse-grained potential
will approximate the many-body potential of mean force. The choice of a
particular method is system-specific and requires a thorough consistency
check. It is important to keep in mind that coarse-graining should be
used with understanding and caution, methods should be crossed-checked
with each other as well as with respect to the reference system.

The package consists of two parts: a C++ kernel and a scripting engine.
The kernel is capable of processing atomistic topologies and
trajectories and offers a flexible framework for reading, manipulating
and analyzing topologies and generated by MD/SD/MC sampling
trajectories. It is modular: new file formats can be integrated without
changing the existing code. Currently, an interface for
GROMACS :raw-latex:`\cite{gromacs4}` topologies and trajectories is
provided. The kernel also includes various coarse-graining tools, for
example calculations of probability distributions of bonded and
non-bonded interactions, correlation and autocorrelation functions, and
updates for the coarse-grained pair potential.

The scripting engine is used to steer the iterative procedures. Here the
analysis tools of the package used for sampling (e.g. GROMACStools) can
be integrated into the coarse-graining workflow, if needed. The
coarse-graining workflow itself is controlled by several Extensible
Markup Language (XML) input files, which contain mapping and other
options required for the workflow control. In what follows, these input
files are described.

Before using the package, do not forget to initalize the variables in
the bash or csh (tcsh)

::

      source <csg-installation>/bin/VOTCARC.bash
      source <csg-installation>/bin/VOTCARC.csh

More details as well as several examples can be found in
ref. :raw-latex:`\cite{Ruehle:2009.a}`. Please cite this paper if you
are using the package. Tutorials can be found on the votcahomepage .

Theoretical background
======================

Mapping
-------

The mapping is an operator that establishes a link between the atomistic
and coarse-grained representations of the system. An atomistic system is
described by specifying the values of the Cartesian coordinates and
momenta

.. math::

   \begin{aligned}
   \bm r^n &=& \{\bm r_1,\dots,\bm r_n\}, \\
   \bm p^n &=& \{\bm p_1,\dots,\bm p_n\}.\end{aligned}

 of the :math:`n` atoms in the system. [1]_ On a coarse-grained level,
the coordinates and momenta are specified by the positions and momenta
of CG sites

.. math::

   \begin{aligned}
   \bm R^N = \{\bm R_1,\dots,\bm R_N\}, \\
   \bm P^N = \{\bm P_1,\dots,\bm P_N\}.\end{aligned}

 Note that capitalized symbols are used for the CG sites while lower
case letters are used for the atomistic system.

The mapping operator :math:`{\bm c}_I` is defined by a matrix for each
bead :math:`I` and links the two descriptions

.. math::

   \begin{aligned}
    {\bm R}_I &=& \sum_{i=1}^{n}c_{Ii}\bm r_i, \\
    {\bm P}_I &=&
    	M_I \dot{{\bm R}}_I =
   	M_I \sum_{i=1}^{n}c_{Ii} \dot{{\bm r}}_i =
   	M_I \sum_{i=1}^{n} \frac{ c_{Ii}} {m_i} {\bm p}_i .
   \label{eq:mapping_scheme}\end{aligned}

 for all :math:`I = 1,\dots,N`.

If an atomistic system is translated by a constant vector, the
corresponding coarse-grained system is also translated by the same
vector. This implies that, for all :math:`I`,

.. math:: \sum_{i=1}^{n}c_{Ii}=1.

In some cases it is useful to define the CG mapping in such a way that
certain atoms belong to several CG beads at the same
time :raw-latex:`\cite{Fritz:2009}`. Following
ref. :raw-latex:`\cite{Noid:2008.1}`, we define two sets of atoms for
each of the :math:`N` CG beads. For each site :math:`I`, a set of
*involved* atoms is defined as

.. math:: {\cal I}_I=\{i|c_{Ii}\ne0\}.

 An atom :math:`i` in the atomistic model is involved in a CG site, *I*,
if and only if this atom provides a nonzero contribution to the sum in
eq. [eq:mapping\_scheme].

A set of *specific* atoms is defined as

.. math:: {\cal S}_I=\{i|c_{Ii}\ne0 \text{ and } c_{Ji}=0 \text{ for all } J \ne I\}.

 In other words, atom :math:`i` is specific to site :math:`I` if and
only if this atom is involved in site :math:`I` and is not involved in
the definition of any other site.

The CG model will generate an equilibrium distribution of momenta that
is consistent with an underlying atomistic model if all the atoms are
*specific* and if the mass of the :math:`I^\text{th}` CG site is given
by :raw-latex:`\cite{Noid:2008.1}`

.. math::

   M_I= \left( \sum_{i \in {\cal I}_I}\frac{c_{Ii}^2}{m_i} \right)^{-1}.
   \label{eq:cg_mass}

 If all atoms are specific and the center of mass of a bead is used for
mapping, then :math:`c_{Ii} = \frac{m_i}{M_I}`, and the
condition [eq:cg\_mass] is automatically satisfied.

Boltzmann inversion
-------------------

Boltzmann inversion is mostly used for *bonded* potentials, such as
bonds, angles, and torsions :raw-latex:`\cite{Tschoep:1998}`. Boltzmann
inversion is structure-based and only requires positions of atoms.

The idea of Boltzmann inversion stems from the fact that in a canonical
ensemble *independent* degrees of freedom :math:`q` obey the Boltzmann
distribution, i. e.

.. math::

   P(q) = Z^{-1} \exp\left[ - \beta U(q) \right]~,
     \label{eq:boltzmann}

 where is a partition function, . Once :math:`P(q)` is known, one can
obtain the coarse-grained potential, which in this case is a potential
of mean force, by inverting the probability distribution :math:`P(q)` of
a variable :math:`q`, which is either a bond length, bond angle, or
torsion angle

.. math::

   U(q) = - k_\text{B} T \ln  P(q) ~.
     \label{eq:inv_boltzmann}

 The normalization factor :math:`Z` is not important since it would only
enter the coarse-grained potential :math:`U(q)` as an irrelevant
additive constant.

Note that the histograms for the bonds :math:`H_r(r)`, angles
:math:`H_\theta(\theta)`, and torsion angles :math:`H_\varphi(\varphi)`
have to be rescaled in order to obtain the volume normalized
distribution functions :math:`P_r(r)`, :math:`P_\theta(\theta)`, and
:math:`P_\varphi(\varphi)`, respectively,

.. math::

   \begin{aligned}
       P_r(r) = \frac{H_r(r)}{4\pi r^2}~,\;
       P_\theta(\theta) = \frac{H_\theta(\theta)}{\sin \theta}~,\;
       P_\varphi(\varphi) = H_\varphi (\varphi)~,
       \label{eq:boltzmann_norm}\end{aligned}

 where :math:`r` is the bond length :math:`r`, :math:`\theta` is the
bond angle, and :math:`\varphi` is the torsion angle. The bonded
coarse-grained potential can then be written as a sum of distribution
functions

.. math::

   \begin{aligned}
       \label{eq:boltzmann_pmf}
       U({r}, \theta, \varphi) &= U_r({r}) + U_{\theta}(\theta) + U_{\varphi}(\varphi)~, \\
       U_q({q}) &= - k_\text{B} T \ln P_q( q ),\; q=r, \theta, \varphi~.
       \nonumber\end{aligned}

On the technical side, the implementation of the Boltzmann inversion
method requires *smoothing* of :math:`U(q)` to provide a continuous
force. Splines can be used for this purpose. Poorly and unsampled
regions, that is regions with high :math:`U(q)`, shall be
*extrapolated*. Since the contribution of these regions to the canonical
density of states is small, the exact shape of the extrapolation is less
important.

Another crucial issue is the cross-correlation of the coarse-grained
degrees of freedom. Independence of the coarse-grained degrees of
freedom is the main assumption that allows factorization of the
probability distribution and the potential, eq. [eq:boltzmann\_pmf].
Hence, one has to carefully check whether this assumption holds in
practice. This can be done by performing coarse-grained simulations and
comparing cross-correlations for all pairs of degrees of freedom in
atomistic and coarse-grained resolution, e. g. using a two-dimensional
histogram, analogous to a Ramachandran plot.  [2]_

Separation of bonded and non-bonded degrees of freedom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When coarse-graining polymeric systems, it is convenient to treat bonded
and non-bonded interactions separately :raw-latex:`\cite{Tschoep:1998}`.
In this case, sampling of the atomistic system shall be performed on a
special system where non-bonded interactions are artificially removed,
so that the non-bonded interactions in the reference system do not
contribute to the bonded interactions of the coarse-grained model.

This can be done by employing exclusion lists using with the option
``—excl``. This is described in detail in sec. [sec:exclusions].

.. figure:: fig/excl
   :alt: [fig:excl]Example of excluded interactions.
   :width: 70.0%

   [fig:excl]Example of excluded interactions.

Iterative methods
-----------------

ht3cm |image|

Iterative workflow control is essential for the IBIand IMCmethods. The
general idea of iterative workflow is sketched in
fig. [fig:iterative\_methods]. A run starts with an initial guess during
the global initialization phase. This guess is used for the first
sampling step, followed by an update of the potential. The update itself
often requires additional postprocessing such as smoothing,
interpolation, extrapolation or fitting. Different methods are available
to update the potential, for instance Iterative Boltzmann Inversion (see
next section [sec:ibi]) or Inverse Monte Carlo (see section [sec:imc]).
The whole procedure is then iterated until a convergence criterion is
satisfied.

Iterative Boltzmann Inversion
-----------------------------

Iterative Boltzmann inversion (IBI) is a natural extension of the
Boltzmann inversion method. Since the goal of the coarse-grained model
is to reproduce the distribution functions of the reference system as
accurately as possible, one can also iteratively refine the
coarse-grained potentials using some numerical scheme.

In IBIthe potential update :math:`\Delta U` is given
by :raw-latex:`\cite{Reith:2003}`

.. math::

   \begin{aligned}
     \label{eq:iter_boltzmann}
     U^{(n+1)} &=& U^{(n)} + \lambda \Delta U^{(n)}~, \\
     \Delta U^{(n)} &=&  k_\text{B} T \ln  \frac{P^{(n)}}{P_{\rm ref}}
     =  U_\text{PMF}^\text{ref} - U_\text{PMF}^{(n)}~.\end{aligned}

 Here :math:`\lambda \in (0,1]` is a numerical factor which helps to
stabilize the scheme.

The convergence is reached as soon as the distribution function
:math:`P^{(n)}` matches the reference distribution function
:math:`P_{\rm ref}`, or, in other words, the potential of mean force,
:math:`U_\text{PMF}^{(n)}`, converges to the reference potential of mean
force.

IBIcan be used to refine both bonded and non-bonded potentials. It is
primarily used for simple fluids with the aim to reproduce the radial
distribution function of the reference system in order to obtain
non-bonded interactions. On the implementation side, IBIhas the same
issues as the inverse Boltzmann method, i. e. smoothing and
extrapolation of the potential must be used.

Inverse Monte Carlo
-------------------

Inverse Monte Carlo (IMC) is an iterative scheme which additionally
includes cross correlations of distributions. A detailed derivation of
the IMCmethod can be found in ref. :raw-latex:`\cite{Lyubartsev:1995}`.

The potential update :math:`\Delta U` of the IMCmethod is calculated by
solving a set of linear equations

.. math::

   \begin{aligned}
       \left<S_{\alpha}\right> - S_{\alpha}^{\text{ref}}= A_{\alpha \gamma} \Delta U_{\gamma}~,
     \label{eq:imc}\end{aligned}

 where

.. math::

   \begin{aligned}
     \label{eq:covariance}
     A_{\alpha \gamma} = \frac{\partial \left< S_{\alpha} \right> }{\partial U_{\gamma}}  =
     \beta \left( \left<S_{\alpha} \right>\left<S_{\gamma} \right> - \left<S_{\alpha} S_{\gamma} \right>  \right)~,
     \nonumber\end{aligned}

 and :math:`S` the histogram of a coarse-grained variable of interest.
For example, in case of coarse-graining of the non-bonded interactions
which depend only on the distance :math:`r_{ij}` between particles
:math:`i` and :math:`j` and assuming that the interaction potential is
short-ranged, i.e. :math:`U(r_{ij})=0` if
:math:`r_{ij} \ge r_{\text{cut} }`, the average value of
:math:`S_{\alpha}` is related to the radial distribution function
:math:`g(r_{\alpha})` by

.. math::

   \left< S_{\alpha} \right> =  \frac{N(N-1)}{2} \frac{4 \pi r_{\alpha}^2 \Delta r} {V}g(r_{\alpha})~,
     \label{eq:s_mean}

 where :math:`N` is the number of atoms in the system
(:math:`\frac{1}{2} N(N-1)` is then the number of all pairs),
:math:`\Delta r` is the grid spacing, :math:`r_{\text{cut}}/M`,
:math:`V` is the total volume of the system. In other words, in this
particular case the physical meaning of :math:`S_{\alpha}` is the number
of particle pairs with interparticle distances
:math:`r_{ij} = r_{\alpha}` which correspond to the tabulated value of
the potential :math:`U_{\alpha}`.

Regularization of Inverse Monte Carlo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To get a well defined cross correlation matrix,
:math:`A_{\alpha \gamma}`, enough sampling is needed. If there is not
enough smapling or the initial potential guess is far from the real
solution of the inverse problem, the algorithm might not converge to a
stable solution. To overcome this instability problem one could
reformulate equation [eq:covariance] by addition of a penalty term. In
this case the potential update is computed as
follows::raw-latex:`\cite{Murtola:2007}`

.. math::

   \label{eq:regularization}
   \Delta U_\gamma = \arg \min \| A_{\alpha \gamma} \Delta U_\gamma - \left(\left<S_{\alpha}\right> - S_{\alpha}^{\text{ref}}\right) \|^2 + \lambda \| R \Delta U_{\gamma} \|^{2}

 Equation [eq:regularization] is known as Tikhonov regularization, where
:math:`R` is the regularization operator, which here is the identity
matrix and :math:`\lambda >0 ` is the regularization parameter. The
optimal choice for :math:`\lambda` can only be determined if the exact
solution of the inverse problem is known, which in practice is not the
case. To get a good initial guess on the magnitude of the regularization
parameter a singular value decomposition of the matrix
:math:`A_{\alpha \gamma}` might help. A good :math:`\lambda` parameter
should dominate the smallest singular values (squared) but is itself
small compared to the larger ones.:raw-latex:`\cite{Rosenberger:2016}`

Force Matching
--------------

Force matching (FM) is another approach to evaluate corse-grained
potentials :raw-latex:`\cite{Ercolessi:1994,Izvekov:2005,Noid:2007}`. In
contrast to the structure-based approaches, its aim is not to reproduce
various distribution functions, but instead to match the multibody
potential of mean force as close as possible with a given set of
coarse-grained interactions.

The method works as follows. We first assume that the coarse-grained
force-field (and hence the forces) depends on :math:`M` parameters
:math:`g_1,...,g_M `. These parameters can be prefactors of analytical
functions, tabulated values of the interaction potentials, or
coefficients of splines used to describe these potentials.

In order to determine these parameters, the reference forces on
coarse-grained beads are calculated by summing up the forces on the
atoms

.. math::

   {{\ensuremath{{\bm F}}}}_I^\text{ref} = \sum_{j \in {\cal S_I}} \frac{d_{Ii}}{c_{Ii}} {{\ensuremath{{\bm f}}}}_j({{\ensuremath{{\bm r}}}^n}),
     \label{eq:force_mapping}

 where the sum is over all atoms of the CG site *I* (see.
sec. [sec:mapping\_operator]). The :math:`d_{Ij}` coefficients can, in
principle, be chosen arbitrarily, provided that the condition
:math:` \sum_{i=1}^{n}d_{Ii}=1` is
satisfied :raw-latex:`\cite{Noid:2008.1}`. If mapping coefficients for
the forces are not provided, it is assumed that :math:`d_{Ij} = c_{Ij}`
(see also sec. [sec:inputfiles]).

By calculating the reference forces for :math:`L` snapshots we can write
down :math:`N \times L` equations

.. math::

   {{\ensuremath{{\bm F}}}}_{Il}^\text{cg}(g_1, \dots ,g_M)={\ensuremath{{\bm F}}}_{il}^\text{ref},\;
     I=1,\dots,N,\; l=1,\dots,L~.
     \label{eq:fmatch1}

 Here :math:`{{\ensuremath{{\bm F}}}}_{Il}^\text{ref}` is the force on
the bead :math:`I` and :math:`{{\ensuremath{{\bm F}}}}_{Il}^\text{cg} `
is the coarse-grained representation of this force. The index :math:`l`
enumerates snapshots picked for coarse-graining. By running the
simulations long enough one can always ensure that
:math:`M < N \times L`. In this case the set of equations [eq:fmatch1]
is overdetermined and can be solved in a least-squares manner.

:math:`{\bm F}_{il}^\text{cg}` is, in principle, a non-linear function
of its parameters :math:`\{g_i\}`. Therefore, it is useful to represent
the coarse-grained force-field in such a way that
equations ([eq:fmatch1]) become linear functions of :math:`\{g_i\}`.
This can be done using splines to describe the functional form of the
forces :raw-latex:`\cite{Izvekov:2005}`. Implementation details are
discussed in ref. :raw-latex:`\cite{Ruehle:2009.a}`.

Note that an adequate sampling of the system requires a large number of
snapshots :math:`L`. Hence, the applicability of the method is often
constrained by the amount of memory available. To remedy the situation,
one can split the trajectory into blocks, find the coarse-grained
potential for each block and then perform averaging over all blocks.

Relative Entropy
----------------

Relative entropy is a method which quantifies the extent of the
configurational phase-space overlap between two molecular
ensembles :raw-latex:`\cite{Wu2005}`. It can be used as a measure of the
discrepancies between various properties of the CG system’s and the
target all-atom (AA) ensemble. It has been shown by Shell
S. :raw-latex:`\cite{Shell2008}` that one can minimize the relative
entropy metric between the model CG system and the target AA system to
optimize CG potential parameters such that the CG ensemble would mimic
the target AA ensemble.

Relative entropy, :math:`S_{\text{rel}}`, is defined as
:raw-latex:`\cite{Shell2008}`

.. math::

   \label{eq:srel}
   S_{\text{rel}} = \sum_{i}p_{\text{AA}}(r_i) \ln\left(
     \frac{p_{\text{AA}}(r_i)}{p_{\text{CG}}\left(M(r_i)\right)}\right) +
   \langle S_{\text{map}} \rangle_{\text{AA}},

 where the sum is over all the configurations of the reference AA
system, :math:`r=\{r_i\} (i=1,2,...)`, :math:`M` is the mapping
operation to generate a corresponding CG configuration, :math:`R_I`,
from a AA configuration, :math:`r_i`, i.e., :math:`R_I = M(r_i)`,
:math:`p_\text{AA}` and :math:`p_\text{CG}` are the configurational
probabilities based on the AA and CG potentials, respectively, and
:math:` \langle
S_{\text{map}}\rangle_{\text{AA}}` is the mapping entropy due to the
average degeneracy of AA configurations mapping to the same CG
configuration, given by

.. math::

   \label{eq:smap}
   S_{\text{map}}(R_I)=\ln\sum_{i}\delta_{R_I,M(r_i)} ,

 where :math:`\delta` is the Kronecker delta function. Physically,
:math:`S_{\text{rel}}` can be interpreted as the likelihood that one
test configuration of the model CG ensemble is representative of the
target AA ensemble, and when the likelihood is a maximum,
:math:`S_{\text{rel}}` is at a minimum. Hence, the numerical
minimization of :math:`S_{\text{rel}}` with respect to the parameters of
the CG model can be used to optimize the CG model.

In a canonical ensemble, substituting canonical configurational
probabilities into eq. [eq:srel], the relative entropy simplifies to

.. math::

   \label{eq:srelcan}
   S_{\text{rel}}=\beta\langle U_{\text{CG}} - U_{\text{AA}}\rangle_{\text{AA}}
   - \beta\left( A_{\text{CG}} - A_{\text{AA}}\right)
   + \langle S_{\text{map}}\rangle_{\text{AA}} ,

 where :math:`\beta={1}/{k_{\text{B}}T}`, :math:`k_{\text{B}}` is the
Boltzmann constant, :math:`T` is the temperature, :math:`U_\text{CG}`
and :math:`U_\text{AA}` are the total potential energies from the CG and
AA potentials, respectively, :math:`A_\text{CG}` and :math:`A_\text{AA}`
are the configurational part of the Helmholtz free energies from the CG
and AA potentials, respectively, and all the averages are computed in
the reference AA ensemble.

Consider a model CG system defined by the CG potentials between various
CG sites such that the CG potentials depend on the parameters
:math:`\boldsymbol\lambda=\{\lambda_1,\lambda_2,...\lambda_n\}`. Then
:math:`\boldsymbol\lambda` are optimized by the relative entropy
minimization. We use the Newton-Raphson strategy for the relative
entropy minimization described in
ref. :raw-latex:`\cite{Chaimovich2011}`. In this strategy, the CG
potential parameters, :math:`\boldsymbol\lambda`, are refined
iteratively as

.. math::

   \label{eq:newtraph}
   \boldsymbol{\lambda} ^{k+1} = \boldsymbol{\lambda} ^{k} -
   \chi \mathbf{H} ^{-1}\cdot
   \nabla_{\lambda} S_{\text{rel}} ,

 where :math:`k` is the iteration index, :math:`\chi\in(0...1)` is the
scaling parameter that can be adjusted to ensure convergence,
:math:`\nabla_{\lambda}S_{\text{rel}}` is the vector of the first
derivatives of :math:`S_{\text{rel}}` with respect to
:math:`\boldsymbol\lambda`, which can be computed from eq. [eq:srelcan]
as

.. math::

   \label{eq:dsrel}
   \nabla_{\lambda}S_{\text{rel}} = \beta \left\langle \frac{\partial
     U_{\text{CG}}}{\partial\lambda}\right\rangle_{\text{AA}} - \beta\left\langle
   \frac{\partial U_{\text{CG}}}{\partial\lambda}\right\rangle_{\text{CG}} ,

 and :math:`\mathbf{H}` is the Hessian matrix of :math:`S_{\text{rel}}`
given by

.. math::

   \begin{aligned}
   \label{eq:Hsrel}
   \mathbf{H}_{ij}&=&\beta \left\langle \frac{\partial^2
     U_{\text{CG}}}{\partial\lambda_i\partial\lambda_j}\right \rangle_{\text{AA}} -
   \beta \left\langle \frac{\partial^2
     U_{\text{CG}}}{\partial\lambda_i\partial\lambda_j}\right \rangle_{\text{CG}}
   \nonumber\\ &&+ \beta^2 \left\langle \frac{\partial
     U_{\text{CG}}}{\partial\lambda_i} \frac{\partial
     U_{\text{CG}}}{\partial\lambda_j}\right\rangle_{\text{CG}} \nonumber\\ &&-
   \beta^2 \left\langle \frac{\partial
     U_{\text{CG}}}{\partial\lambda_i}\right\rangle_{\text{CG}} \left\langle
   \frac{\partial U_{\text{CG}}}{\partial\lambda_j}\right\rangle_{\text{CG}}.\end{aligned}

To compute :math:`\nabla_{\lambda}S_{\text{rel}}` and :math:`\mathbf{H}`
from eq. [eq:dsrel] and [eq:Hsrel], we need average CG energy
derivatives in the AA and CG ensembles. For two-body CG pair potentials,
:math:`u_{\text{CG}}`, between CG sites, the ensemble averages of the CG
energy derivatives can be computed as

.. math::

   \begin{aligned}
   \left\langle\left(\frac{\partial^a U_{\text{CG}}}{\partial \lambda^a}\right)^b
   \right\rangle_{\text{AA}}& =
   &\left\langle\left(\sum_{i<j}\frac{\partial^{a}u_{\text{CG}}(r_{ij})}
   {\partial \lambda^a}\right)^b\right\rangle_{\text{AA}}\nonumber \\
   \left\langle\left(\frac{\partial^a U_{\text{CG}}}{\partial \lambda^a}\right)^b
   \right\rangle_{\text{CG}}& =
   &\left\langle\left(\sum_{i<j}\frac{\partial^{a}u_{\text{CG}}(r_{ij})}
   {\partial \lambda^a}\right)^b\right\rangle_{\text{CG}}  ,\end{aligned}

 where the sum is performed over all the CG site pairs :math:`(i,j)`,
:math:`a` stands for the 1\ :math:`^{\text{st}}`,
2\ :math:`^{\text{nd}}`,... derivatives and :math:`b` stands for the
different powers, i.e., :math:`b=1,2,...`. For the averages in the AA
ensemble, first a single AA system simulation can be performed and RDFs
between the CG sites in the AA ensemble can be saved, then the average
CG energy derivatives in AA ensemble can be computed by processing the
CG RDFs in the AA ensemble using the CG potentials at each iteration.
For the averages in the CG ensemble, since the CG ensemble changes with
the CG parameters, :math:`\boldsymbol\lambda`, a short CG simulation is
performed at each iteration to generate corresponding CG configurations.

Comparisons between relative entropy and other coarse-graining methods
are made in ref. :raw-latex:`\cite{rudzinski_coarse-graining_2011}`
and :raw-latex:`\cite{Chaimovich2011}`. Chaimovich and
Shell :raw-latex:`\cite{Chaimovich2011}` have shown that for certain CG
models relative entropy minimization produces the same CG potentials as
other methods, e.g., it is equivalent to the IBI when CG interactions
are modeled using finely tabulated pair additive potentials, and to the
FM when a CG model is based on :math:`N-`\ body interactions, where
:math:`N` is the number of degrees of freedom in the CG model. However,
there are some advantages of using relative entropy based
coarse-graining. Relative entropy method allows to use analytical
function forms for CG potentials, which are desired in theoretical
treatments, such as parametric study of CG potentials, whereas, methods,
like IBI, use tabulated potentials. Recently Lyubartsev et.
al :raw-latex:`\cite{lyubartsev2010systematic}` have shows how to use
IMC with an analytical function form, too. BI, IBI, and IMC methods are
based on pair correlations and hence, they are only useful to optimize
2-body CG potentials, whereas, relative entropy uses more generic metric
which offers more flexibility in modeling CG interactions and not only
2-body, but also 3-body (for example see
ref. :raw-latex:`\cite{lu_coarse-graining_2014}`) and N-body CG
potentials can be optimized. In addition to the CG potential
optimization, the relative entropy metric can also be used to optimize
an AA to CG mapping operator.

Input files
===========

Mapping files
-------------

Mapping relates atomistic and coarse-grained representations of the
system. It is organized as follows: for each molecule *type* a mapping
file is created. When used as a command option, these files are combined
in a list separated by a semicolon, e. g.
``—cg`` ``protein.xml;solvent.xml``.

ht6cm |image|

Each mapping file contains a *topology* of the coarse-grained molecule
and a list of *maps*. Topology specifies coarse-grained beads and bonded
interactions between them. Each coarse-grained bead has a name, type, a
list of atoms which belong it, and a link to a map. A map is a
:math:`c_{Ii}` for an atom :math:`i` belonging to the bead :math:`I`. It
is used to calculate the position of a coarse-grained bead from the
positions of atoms which belong to it. Note that :math:`c_{Ii}` will be
automatically re-normalized if their sum is not equal to 1, i. e. in the
case of a center-of-mass mapping one can simply specify atomic masses. A
complete reference for mapping file definitions can be found in
sec. [sec:ref\_mapping].

As an example, we will describe here a mapping file of a united atom
model of a propane molecule, chemical structure of which is shown in
fig. [fig:intro:propane]. In this coarse-grained model two bead types
(A,B) and three beads (A1, B1, A2) are defined, as shown in
fig. [fig:propane\_map]. We will use centers of mass of the beads as
coarse-grained coordinates.

Extracts from the ``propane.xml`` file of the tutorial are shown below.
The tag indicates the molecule name in the coarse-grained topology. The
tag must match the name of the molecule in the atomistic representation.
In the section all beads are defined by specifying bead name (A1, B1,
A2), type, and atoms belonging to this bead in the form
``residue id:residue name:atom name``. For each bead a map has to be
specified, which is defined later in section. Note that bead and can be
different, which might be useful in a situation when chemically
different beads (A1, B1) are assigned to the same bead type. After
defining all beads the bonded interactions of the coarse-grained
molecule must be specified in the section. This is done by using the
identifiers of the beads in the coarse-grained model. Finally, in the
section, the mapping coefficients are defined. This includes a weighting
of the atoms in the topology section. In particular, the number of
weights given should match the number of beads.

Verification of a mapping
-------------------------

Note that the tag should match the molecule name in the reference
system. A common mistake is that beads have wrong names. In this case,
the tool can be used in order to identify the atoms which are read in
from a topology file ``.tpr``. This tool displays the atoms in the
format ``residue id:residue name:atom name``. For multicomponent
systems, it might happen that molecules are not identified correctly.
The workaround for this case is described in sec. [sec:adv\_topology].

To compare coarse-grained and atomistic configurations one can use a
standard visualization program, e. g. ``vmd``. When comparing
trajectories, one has to be careful, since ``vmd`` opens both a ``.gro``
and ``.trr`` file. The first frame is then the ``.gro`` file and the
rest is taken from the ``.trr`` file. The coarse-grained trajectory
contains only the frames of the trajectory. Hence, the first frame of
the atomistic run has to be removed using the ``vmd`` menu.

Advanced topology handling
--------------------------

A topology is completely specified by a set of beads, their types, and a
list of bonded interactions. votcais able to read topologies in the
GROMACS\ ``.tpr`` format. For example, one can create a coarse-grained
topology based on the mapping file and atomistic GROMACStopology using .

::

      csg_gmxtopol --top topol.tpr --cg propane.xml --out out.top

In some cases, however, one might want to use a .pdb, H5MD or .dump file
which does not contain all information about the atomistic topology. In
this case, additional information can be supplied in the XMLmapping
file.

A typical example is lack of a clear definition of molecules, which can
be a problem for simulations with several molecules with multiple types.
During coarse-graining, the molecule type is identified by a name tag as
names must be clearly identified. To do this, it is possible to read a
topology and then modify parts of it. The new XMLtopology can be used
with the ``—tpr`` option, as any other topology file.

For example, if information about a molecule is not present at all, one
can create one from a ``.pdb`` file as follows

::

    <topology base="snapshot.pdb">
      <molecules>
        <clear/>
        <define name="mCP" first="1" nbeads="52" nmols="216"/>
      </molecules>
    </topology>

where :math:`<`\ clear/\ :math:`>` clears all information that was
present before.

Old versions of GROMACSdid not store molecule names. In order to use
this feature, a recent ``.tpr`` file containing molecule names should
always be provided. For old topologies, rerun GROMACS to update the old
topology file.

If molecule information is already present in the parent topology but
molecules are not named properly (as it is the case with old
GROMACS\ ``.tpr`` files), one can rename them using

::

     <topology base="topol.tpr">
       <molecules>
         <rename name="PPY3" range="1:125"/>
         <rename name="Cl" range="126:250"/>
       </molecules>
     </topology>

Here, the file ``topol.tpr`` is loaded first and all molecules are
renamed afterwards.

If you do not have a .pdb/.gro file and you want to read trajectory from
LAMMPS .dump file or H5MD then it is also possible to directly define
topology in XMLfile. Here is an example of such file where the
trajectory is read from H5MD file:

::

    <topology>
        <!-- particle group name in H5MD file -->
        <h5md_particle_group name="atoms" />
        <molecules>
            <!-- define molecule, number of beads, number of mols -->
            <molecule name="BUT" nmols="4000" nbeads="4">
                <!-- composition of molecule, bead definition -->
                <bead name="C1" type="C" mass="15.035" q="0.0" />
                <bead name="C2" type="C" mass="14.028" q="0.0" />
                <bead name="C3" type="C" mass="14.028" q="0.0" />
                <bead name="C4" type="C" mass="15.035" q="0.0" />
            </molecule>
        </molecules>
        <!-- bonded terms -->
        <bonded>
            <bond>
                <name>bond1</name>
                <beads>
                    BUT:C1 BUT:C2
                </beads>
            </bond>
            <bond>
                <name>bond2</name>
                <beads>
                    BUT:C2 BUT:C3
                </beads>
            </bond>
            <angle>
                <name>angle1</name>
                <beads>
                    BUT:C1 BUT:C2 BUT:C3
                    BUT:C2 BUT:C3 BUT:C4
                </beads>
            </angle>
            <dihedral>
                <name>dihedral1</name>
                <beads>
                    BUT:C1 BUT:C2 BUT:C3 BUT:C4
                </beads>
            </dihedral>
        </bonded>
    </topology>

The list of molecules is defined in section ``molecules`` where every
``molecule`` is replicated ``nmols`` times. Inside ``molecule`` the list
of ``bead`` has to be defined with the name, type, mass and charge.

The box size can be set by the tag ``box``:

::

    <box xx="6.0" yy="6.0" zz="6.0" />

where ``xx, yy, zz`` are the dimensions of the box.

A complete reference for XMLtopology file can be found in
sec. [sec:ref\_topology].

Trajectories
------------

A trajectory is a set of frames containing coordinates (velocities and
forces) for the beads defined in the topology. votcacurrently supports
``.trr``, ``.xtc``, ``.pdb``, ``.gro`` and H5MD ``.h5`` trajectory
formats.

Once the mapping file is created, it is easy to convert an atomistic to
a coarse-grained trajectory using

::

      csg_map --top topol.tpr --trj traj.trr --cg propane.xml --out cg.gro

The program also provides the option ``—no-map``. In this case, no
mapping is done and works as a trajectory converter. In general, mapping
can be enabled and disabled in most analysis tools, e.g. in or .

Note that the topology files can have a different contents as bonded
interactions are not provided in all formats. In this case, mapping
files can be used to define and relabel bonds.

Also note that the default setting concerning mapping varies
individually between the programs. Some have a default setting that does
mapping (such as , use ``—no-map`` to disable mapping) and some have
mapping disabled by default (e.g. , use ``—cg`` to enable mapping).

Setting files
-------------

::

    <cg>
      <non-bonded> <!-- non-bonded interactions -->
        <name>A-A</name> <!-- name of the interaction -->
        <type1>A</type1> <!-- types involved in this interaction -->
        <type2>A</type2>
        <min>0</min>  <!-- dimension + grid spacing of tables-->
        <max>1.36</max>
        <step>0.01</step>
        <inverse>
          ... specific commands
        </inverse>

        ... specific section for inverse boltzmann, force matching etc.
      </non-bonded>
    </cg>

A setting file is written in the format ``.xml``. It consists of a
general section displayed above, and a specific section depending on the
program used for simulations. The setting displayed above is later
extended in the sections on iterative boltzmann inversion (), force
matching () or statistical analysis ().

Generally, is an analysis tool which can be used for computing radial
distribution functions and analysing them. As an example, the command

::

      csg_stat --top topol.tpr --trj traj.xtc --options settings.xml

computes the distributions of all interactions specified in
``settings.xml`` and writes all tabulated distributions as files
``interaction name.dist.new``.

Table formats
-------------

In the iterative framework distribution functions, potentials and forces
are returned as tables and saved in a file. Those tables generally have
the format

::

      x y [error] flag

where ``x`` is input quantity (e.g. radius :math:`r`, angles
:math:`\theta` or :math:`\phi`), ``y`` is the computed quantity (e.g. a
potential) and ``[error]`` is an optional error for ``y``. The token
``flag`` can take the values ``i``, ``o`` or ``u``. In the first case,
``i`` (``in range``) describes a value that lies within the data range,
``o`` (``out of range``) symbolises a value out of the data range and
``u`` stands for an ``undefined`` value.

The token ``flag`` will be important when extrapolating the table as
described in sec. [sec:post\_processing].

For historical reasons, uses a slightly different table format, it has
no ``flag`` column and uses the third column as a force column when
outputting a potential.

Preparing coarse-grained runs
=============================

Preliminary note
~~~~~~~~~~~~~~~~

The coarse-grained run requires the molecule topology on the one hand
and suitable potentials on the other. In this chapter, the generation of
coarse-grained runs is decribed next, followed by a post-processing of
the potential.

If the potential is of such a form that it allows direct fitting of a
functional form, the section on post-processing can be skipped. Instead,
a program of choice should be used to fit a functional form to the
potential. Nevertheless, special attention should be paid to units
(angles, bondlengths). The resulting curve can then be specified in the
MD package used for simulation. However, most potentials don’t allow an
easy processing of this kind and tabulated potentials have to be used.

Generating a topology file for a coarse-grained run
---------------------------------------------------

WARNING: This section describes experimental features. The exact names and options of the program might change in the near future. The section is specific to GROMACSsupport though a generalization for other MD packages is planned.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mapping definition is close to a topology needed for a coarse
grained run. To avoid redundant work, can be used to automatically
generate a gromacs topology based on an atomistic reference system and a
mapping file.

At the current state, can only generate the topology for the first
molecule in the system. If more molecule types are present, a special
tpr file has to be prepared. The program can be executed by

::

      csg_gmxtopol --top topol.tpr --cg map.xml --out cgtop

which will create a file ``cgtop.top``. This file includes the topology
for the first molecule including definitions for atoms, bonds, angles
and dihedrals. It can directly be used as a topology in GROMACS, however
the force field definitions (atom types, bond types, etc.) still have to
be added manually.

Post-processing of the potential
--------------------------------

The votcapackage provides a collection of scripts to handle potentials.
They can be modified, refined, integrated or inter- and extrapolated.
These scripts are the same ones as those used for iterative methods in
chapter [sec:iterative\_methods]. Scripts are called by . A complete
list of available scripts can be found in sec. [sec:csg\_table].

The post-processing roughly consists of the following steps (see further
explanations below):

-  (manually) clipping poorly sampled (border) regions

-  resampling the potential in order to change the grid to the proper
   format ()

-  extrapolation of the potential at the borders ( table extrapolate)

-  exporting the table to xvg ( convert\_potential gromacs)

Clipping of poorly sampled regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Regions with an irregular distribution of samples should be deleted
first. This is simply done by editing the ``.pot`` file and by deleting
those values.

Alternatively, manually check the range where the potential still looks
good and is not to noisy and set the flags in the potential file of the
bad parts by hand to ``o`` (for ``out of range``). Those values will
later be extrapolated and overwritten.

Resampling
~~~~~~~~~~

Use the command

::

      csg_resample --in table.pot --out table_resample.pot \
                   --grid min:step:max

to resample the potential given in file –``table.pot`` from ``min`` to
``max`` with a grid spacing of ``step`` steps. The result is written to
the file specified by ``out``. Additionally, allows the specification of
spline interpolation (``spfit``), the calculation of derivatives
(``derivative``) and comments (``comment``). Check the help (``help``)
for further information.

It is important to note that the values ``min`` and ``max`` *don’t*
correspond to the minimum and maximum value in the input file, but to
the range of values the potential is desired to cover after
extrapolation. Therefore, values in :math:`[ \min,\max ]` that are not
covered in the file are automatically marked by a flag ``o`` (for
``out of range``) for extrapolation in the next step.

The potential *don’t* have to start at 0, this is done by the export
script (to xvg) automatically.

Extrapolation
~~~~~~~~~~~~~

The following line

::

      csg_call table extrapolate [options] table_resample.pot \
               table_extrapolate.pot

calls the extrapolation procedure, which processes the range of values
marked by . The input file is ``table_resample.pot`` created in the last
step.

After resampling, all values in the potential file that should be used
as a basis for extrapolation are marked with an ``i``, while all values
that need extrapolation are marked by ``o``. The command above now
extrapolates all ``o`` values from the ``i`` values in the file.
Available options include averaging over a certain number of points
(``avgpoints``), changing the functional form (``function``, default is
quadratic), extrapolating just the left or right region of the file
(``region``) and setting the curvature (``curvature``).

The output ``table_extrapolate.pot`` of the extrapolation step can now
be used for the coarse-grained run. If GROMACSis used as a molecule
dynamics package, the potential has to be converted and exported to a
suitable GROMACSformat as described in the final step.

Exporting the table
~~~~~~~~~~~~~~~~~~~

Finally, the table is exported to ``xvg``. The conversion procedure
requires a small xml file ``table.xml`` as shown below:

::

      <cg>
        <non-bonded>
          <name>XXX</name>
          <step>0.01</step>
        </non-bonded>
        <inverse>
          <gromacs>
            <pot_max>1e8</pot_max>
            <table_end>8.0</table_end>
            <table_bins>0.002</table_bins>
          </gromacs>
        </inverse>
      </cg>

where ``<table_end>`` is the GROMACS\ ``rvdw+table_extension`` and
``<pot_max>`` is just a number slightly smaller than the upper value of
single/ double precision. The value given in ``<table_bins>``
corresponds to the ``step`` value of
``csg_resample -grid min:step:max``.

Using the ``xml`` file above, call

::

      csg_call --options table.xml --ia-type non-bonded --ia-name XXX \
        convert_potential gromacs table_extrapolate.pot table.xvg

to convert the extrapolated values in ``table_extrapolate.pot`` to
``table.xvg`` (The file will contain the GROMACSC12 parts only which are
stored in the sixth und seventh column, this can be changed by adding
the ``–ia-type C6`` option (for the fourth and fiveth column) or
``–ia-type CB`` option (for the second and third column) after . Ensure
compatibility with the GROMACStopology. See the GROMACSmanual for
further information).

To obtain a bond table, run

::

      csg_call --ia-type bond --ia-name XXX --options table.xml \
      convert_potential gromacs table_extrapolate.pot table.xvg

It is also possible to use ``angle`` and ``dihedral`` as type as well,
but make to sure to have a ``bonded`` section similar to the
``non-bonded`` section above with the corresponding interaction name.

Internally ``convert_potential gromacs`` will do the following steps:

-  Resampling of the potential from 0 (or -180 for dihedrals) to
   ``table_end`` (or 180 for angles and dihedrals) with step size
   ``table_bins``. This is needed for gromacs the table must start with
   0 or -180.

-  Extrapolate the left side (to 0 or -180) exponentially

-  Extrapolate the right side (to ``table_end`` or 180) exponentially
   (or constant for non-bonded interactions)

-  Shift it so that the potential is zero at ``table_end`` for
   non-bonded interactions or zero at the minimum for bonded interaction

-  Calculate the force (assume periodicity for dihedral potentials)

-  Write to the format needed by gromacs

An example on non-bonded interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

      csg_call pot shift_nonbonded table.pot table.pot.refined
      csg_resample --grid 0.3:0.05:2 --in table.pot.refined \
               --out table.pot.refined
      csg_call table extrapolate --function quadratic --region left \
               table.pot.refined table.pot.refined
      csg_call table extrapolate --function constant --region right \
               table.pot.refined table.pot.refined

Alternatives
------------

Additionally to the two methods described above, namely (a) providing
the MD package directly with a functional form fitted with a program of
choice or (b) using ``csg_resample``, ``csg_call table extrapolate`` and
``csg_call convert_potential``, another method would be suitable. This
is integrating the force table as follows

::

      -Integrate the table
      $csg_call table integrate force.d minus_pot.d
      -multiply by -1
      $csg_call table linearop minus_pot.d pot.d -1 0

Boltzmann Inversion
===================

provides a potential of mean force for a given degree of freedom.

t7cm |image|

It is mostly used for deriving *bonded* interactions from canonical
sampling of a single molecule in vacuum, e. g. for polymer
coarse-graining, where it is difficult to separate bonded and non-bonded
degrees of freedom :raw-latex:`\cite{Tschoep:1998}`. The non-bonded
potentials can then be obtained by using iterative methods or force
matching.

The main tool which can be used to calculate histograms, cross-correlate
coarse-grained variables, create exclusion lists, as well as prepare
tabulated potentials for coarse-grained simulations is . It parses the
whole trajectory and stores all information on bonded interactions in
memory, which is useful for interactive analysis. For big systems,
however, one can run out of memory. In this case can be used which,
however, has a limited number of tasks it can perform (see sec.
[sec:setting\_files] for an example on its usage).

Another useful tool is . It can be used to convert an atomistic
trajectory to a coarse-grained one, as it is discussed in
sec. [sec:trajectory].

To use one has to first define a mapping scheme. This is outlined
in sec. [sec:mapping\_files]. Once the mapping scheme is specified, it
is possible to generate an exclusion list for the proper sampling of the
atomistic resolution system.

Generating exclusion lists
--------------------------

Exclusion lists are useful when sampling from a special reference system
is needed, for example for polymer coarse-graining with a separation of
bonded and non-bonded degrees of freedom.

To generate an exclusion list, an atomistic topology without exclusions
and a mapping scheme have to be prepared first. Once the .tpr topology
and .xml mapping files are ready, simply run

::

      csg_boltzmann --top topol.tpr --cg mapping.xml --excl exclusions.txt

This will create a list of exclusions for all interactions that are not
within a bonded interaction of the coarse-grained sub-bead. As an
example, consider coarse-graining of a linear chain of three beads which
are only connected by bonds. In this case, will create exclusions for
all non-bonded interactions of atoms in the first bead with atoms of the
3rd bead as these would contribute only to the non-bonded interaction
potential. Note that will only create the exclusion list for the fist
molecule in the topology.

To add the exclusions to the GROMACStopology of the molecule, either
include the file specified by the –excl option into the .top file as
follows

::

      [ exclusions ]
      #include "exclusions.txt"

or copy and paste the content of that file to the exclusions section of
the gromacs topology file.

Statistical analysis
--------------------

For statistical analysis provides an interactive mode. To enter the
interactive mode, use the ``–trj`` option followed by the file name of
the reference trajectory

::

      csg_boltzmann --top topol.tpr --trj traj.trr --cg mapping.xml

To get help on a specific command of the interactive mode, type

::

      help <command>

for example

::

      help hist
      help hist set periodic

Additionally, use the

::

      list

command for a list of available interactions. Note again that loads the
whole trajectory and all information on bonded interactions into the
memory. Hence, its main application should be single molecules. See the
introduction of this chapter for the command.

If a specific interaction shall be used, it can be referred to by

::

      molecule:interaction-group:index

Here, ``molecule`` is the molecule number in the whole topology,
``interaction-group`` is the name specified in the ``<bond>`` section of
the mapping file, and ``index`` is the entry in the list of
interactions. For example, ``1:AA-bond:10`` refers to the 10th bond
named ``AA-bond`` in molecule 1. To specify a couple of interactions
during analysis, either give the interactions separated by a space or
use wildcards (e.g. ``*:AA-bond*``).

To exit the interactive mode, use the command ``q``.

If analysis commands are to be read from a file, use the pipe or stdin
redirects from the shell.

::

      cat commands | csg_boltzmann topol.top --trj traj.trr --cg mapping.xml

Distribution functions and tabulated potentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Distribution functions (tabulated potentials) can be created with the
``hist`` (``tab``) command. For instance, to write out the distribution
function for all interactions of group AA-bond (where AA-bond is the
name specified in the mapping scheme) to the file AA.txt, type

::

      hist AA.txt *:AA-bond:*

The command

::

      hist set

prints a list of all parameters that can be changed for the histogram:
the number ``n`` of bins for the table, bounds ``min`` and ``max`` for
table values, scaling and normalizing, a flag ``periodic`` to ensure
periodic values in the table and an ``auto`` flag. If ``auto`` is set to
1, bounds are calculated automatically, otherwise they can be specified
by ``min`` and ``max``. Larger values in the table might extend those
bounds, specified by parameter ``extend``.

To directly write the Boltzmann-inverted potential, the ``tab`` command
can be used. Its usage and options are very similar to the ``hist``
command. If tabulated potentials are written, special care should be
taken to the parameters ``T`` (temperature) and the ``scale``. The
``scale`` enables volume normalization as given in
eq. [eq:boltzmann\_norm]. Possible values are ``no`` (no scaling),
``bond`` (normalize bonds) and ``angle`` (normalize angles). To write
out the tabulated potential for an angle potential at a temperature of
300K, for instance, type:

::

      tab set T 300
      tab set scale angle
      tab angle.pot *:angle:*

The table is then written into the file ``angle.pot`` in the format
described in sec. [sec:table\_formats]. An optional correlation analysis
is described in the next section. After the file has been created by
command ``tab``, the potential is prepared for the coarse-grained run in
chapter [sec:usage:cgrun].

Correlation analysis
~~~~~~~~~~~~~~~~~~~~

The factorization of :math:`P` in eq. [eq:boltzmann\_pmf] assumed
uncorrelated quantities. offers two ways to evaluate correlations of
interactions. One option is to use the linear correlation coefficient
(command ``cor``).

However, this is not a good measure since ``cor`` calculates the linear
correlation only which might often lead to misleading
results :raw-latex:`\cite{Ruehle:2009.a}`. An example for such a case
are the two correlated random variables :math:`X \sim U[-1,1]` with
uniform distribution, and :math:`Y:=X^2`. A simple calculation shows
:math:`cov(X,Y)=0` and therefore

.. math:: cor=\frac{cov(X,Y)}{\sqrt{var(X)var(Y)}}=0.

 A better way is to create 2D histograms. This can be done by specifying
all values (e.g. bond length, angle, dihedral value) using the command
*vals*, e.g.:

::

      vals vals.txt 1:AA-bond:1 1:AAA-angle:A

This will create a file which contains 3 columns, the first being the
time, and the second and third being bond and angle, respectively.
Columns 2 and 3 can either be used to generate the 2D histogram, or a
simpler plot of column 3 over 2, whose density of points reflect the
probability.

Two examples for 2D histograms are shown below: one for the propane
molecule and one for hexane.

.. figure:: fig/propane_hist2d
   :alt: hexane histograms: before and after the coarse-grained run
   :height: 2.80000cm

   hexane histograms: before and after the coarse-grained run

.. figure:: fig/hexane2
   :alt: hexane histograms: before and after the coarse-grained run
   :height: 3.60000cm

   hexane histograms: before and after the coarse-grained run

The two plots show the correlations between angle and bondlength for
both molecules. In the case of propane, the two quantities are not
correlated as shown by the centered distribution, while correlations
exist in the case of hexane. Moreover, it is visible from the hexane
plot that the partition of the correlations has changed slightly during
coarse-graining.

The tabulated potentials created in this section can be further modified
and prepared for the coarse-grained run: This includes fitting of a
smooth functional form, extrapolation and clipping of poorly sampled
regions. Further processing of the potential is decribed in chapter
[sec:usage:cgrun].

Force matching
==============

.. figure:: usage/fig/flow_fmatch.eps
   :alt: Flowchart to perform force matching.
   :width: 50.0%

   Flowchart to perform force matching.

The force matching algorithm with cubic spline basis is implemented in
the utility. A list of available options can be found in the reference
section of (command ``–h``).

Program input
-------------

needs an atomistic reference run to perform coarse-graining. Therefore,
the trajectory file *must contain forces* (note that there is a suitable
option in the GROMACS\ ``.mdp`` file), otherwise will not be able to
run.

In addition, a mapping scheme has to be created, which defines the
coarse-grained model (see sec. [sec:inputfiles]). At last, a control
file has to be created, which contains all the information for
coarse-graining the interactions and parameters for the force-matching
run. This file is specified by the tag ``–options`` in the XMLformat. An
example might look like the following

::

       <cg>
         <!--fmatch section -->
         <fmatch>
           <!--Number of frames for block averaging -->
           <frames_per_block>6</frames_per_block>
           <!--Constrained least squares?-->
           <constrainedLS>false</constrainedLS>
         </fmatch>
         <!-- example for a non-bonded interaction entry -->
         <non-bonded>
           <!-- name of the interaction -->
           <name>CG-CG</name>
           <type1>A</type1>
           <type2>A</type2>
           <!-- fmatch specific stuff -->
           <fmatch>
             <min>0.27</min>
             <max>1.2</max>
             <step>0.02</step>
             <out_step>0.005</out_step>
           </fmatch>
         </non-bonded>
       </cg>

Similarly to the case of spline fitting (see sec. [sec:ref\_programs] on
), the parameters ``min`` and ``max`` have to be chosen in such a way as
to avoid empty bins within the grid. Determining ``min`` and ``max`` by
using is recommended (see sec. [sec:setting\_files]). A full description
of all available options can be found in sec. [sec:ref\_options].

Program output
--------------

produces a separate ``.force`` file for each interaction, specified in
the CG-options file (option ``options``). These files have 4 columns
containing distance, corresponding force, a table flag and the force
error, which is estimated via a block-averaging procedure. If you are
working with an angle, then the first column will contain the
corresponding angle in radians.

To get table-files for GROMACS, integrate the forces in order to get
potentials and do extrapolation and potentially smoothing afterwards.

Output files are not only produced at the end of the program execution,
but also after every successful processing of each block. The user is
free to have a look at the output files and decide to stop , provided
the force error is small enough.

Integration and extrapolation of .force files 
----------------------------------------------

To convert forces (``.force``) to potentials (``.pot``), tables have to
be integrated. To use the built-in integration command from the
scripting framework, execute

::

     $csg_call table integrate CG-CG.force minus_CG-CG.pot
     $csg_call table linearop minus_CG-CG.d CG-CG.d -1 0

This command calls the script, which integrates the force and writes the
potential to the ``.pot`` file.

In general, each potential contains regions which are not sampled. In
this case or in the case of further post-processing, the potential can
be refined by employing resampling or extrapolating methods. See sec.
[sec:post\_processing] for further details.

Iterative methods
=================

The following sections deal with the methods of Iterative Boltzmann
Inversion (IBI), Inverse Monte Carlo (IMC), and Relative Entropy (RE).

In general, IBI, IMC, and REare implemented within the same framework.
Therefore, most settings and parameters of those methods are similar and
thus described in a general section (see sec.
[sec:iterative\_methods\_imc]). Further information on iterative methods
follows in the next chapters, in particular on the IBI, IMC, and
REmethods.

.. figure:: usage/fig/flow_ibi.eps
   :alt: [fig:flow\_ibi]Flowchart to perform iterative Boltzmann
   inversion.
   :width: 7.00000cm

   [fig:flow\_ibi]Flowchart to perform iterative Boltzmann inversion.

Iterative workflow control
--------------------------

.. figure:: functionality/fig/flowchart.eps
   :alt: [fig:flowchart]Block-scheme of the workflow control for the
   iterative methods. The most time-consuming parts are marked in red.
   :width: 7.00000cm

   [fig:flowchart]Block-scheme of the workflow control for the iterative
   methods. The most time-consuming parts are marked in red.

Iterative workflow control is essential for the IBI, IMC, and REmethods.

The general idea of iterative workflow is sketched in
fig. [fig:flowchart]. During the global initialization the initial guess
for the coarse-grained potential is calculated from the reference
function or converted from a given potential guess into the internal
format. The actual iterative step starts with an iteration
initialization. It searches for possible checkpoints and copies and
converts files from the previous step and the base directory. Then, the
simulation run is prepared by converting potentials into the format
required by the external sampling program and the actual sampling is
performed.

After sampling the phasespace, the potential update is calculated.
Often, the update requires postprocessing, such as smoothing,
interpolation, extrapolation or fitting to an analytical form.

Finally, the new potential is determined and postprocessed. If the
iterative process continues, the next iterative step will start to
initialize.

How to start:
^^^^^^^^^^^^^

The first thing to do is generate reference distribution functions.
These might come from experiments or from atomistic simulations. To get
reasonable results out of the iterative process, the reference
distributions should be of good quality (little noise, etc).

 votcacan create initial guesses for the coarse-grained potentials by
boltzmann inverting the distribution function. If a custom initial guess
for an interaction shall be used instead, the table can be provided in
*:math:`<`\ interaction\ :math:`>`.pot.in*. As already mentioned,
votcaautomatically creates potential tables to run a simulation.
However, it does not know how to run a coarse-grained simulation.
Therefore, all files needed to run a coarse-grained simulation, except
for the potentials that are iteratively refined, must be provided and
added to the in the settings XML-file. If an atomistic topology and a
mapping definition are present, votcaoffers tools to assist the setup of
a coarse-grained topology (see chapter [sec:usage:cgrun]).

To get an overview of how input files look like, it is suggested to take
a look at one of the tutorials provided on .

In what follows we describe how to set up the iterative coarse-graining,
run the main script, continue the run, and add customized scripts.

Preparing the run
~~~~~~~~~~~~~~~~~

To start the first iteration, one has to prepare the input for the
sampling program. This means that all files for running a coarse-grained
simulation must be present and described in a separate XMLfile, in our
case ``settings.xml`` (see sec. [sec:setting\_files] for details). An
extract from this file is given below. The only exception are tabulated
potentials, which will be created and updated by the script in the
course of the iterative process.

The input files include: target distributions, initial guess (optional)
and a list of interactions to be iteratively refined. As a target
distribution, any table file can be given (e.g. GROMACSoutput from
``g_rdf``). The program automatically takes care to resample the table
to the correct grid spacing according to the options provided in
``settings.xml``.

The initial guess is normally taken as a potential of mean force and is
generated by Boltzmann-inversion of the corresponding distribution
function. It is written in ``step_000/<name>.pot.new``. If you want to
manually specify the initial guess for a specific interaction, write the
potential table to a file called ``<name>.pot.in`` in the folder where
you plan to run the iterative procedure.

A list of interactions to be iteratively refined has to be given in the
options file. As an example, the ``setting.xml`` file for a propane is
shown in listing [list:settings]. For more details, see the full
description of all options in ref. [sec:ref\_options].

Starting the iterative process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After all input files have been set up, the run can be started by

::

      csg_inverse --options settings.xml

Each iteration is stored in a separate directory, named
``step_<iteration>``. ``step_000`` is a special folder which contains
the initial setup. For each new iteration, the files required to run the
CG simulation (as specified in the config file) are copied to the
current working directory. The updated potentials are copied from the
last step, ``step_<n-1>/<interaction>.pot.new``, and used as the new
working potentials ``step_<n>/<interaction>.pot.cur``.

After the run preparation, all potentials are converted into the format
of the sampling program and the simulation starts. Once the sampling has
finished, analysis programs generate new distributions, which are stored
in ``<interaction>.dist.new``, and new potential updates, stored in
``<interaction>.dpot.new``.

Before adding the update to the old potential, it can be processed in
the ``post_update`` step. For each script that is specified in the
postupdate, ``<interaction>.dpot.new`` is renamed to
``<interaction>.dpot.old`` and stored in
``<interaction>.dpot.<a-number>`` before the processing script is
called. Each processing script uses the current potential update
``<interaction>.dpot.cur`` and writes the processed update to
``<interaction>.dpot.new``. As an example, a pressure correction is
implemented as a postupdate script within this framework.

After all postupdate scripts have been called, the update is added to
the potential and the new potential ``<interaction>.pot.new`` is
written. Additional post-processing of the potential can be performed in
the ``post_add`` step which is analogous to the ``post_update`` step
except for a potential instead of an update.

To summarize, we list all standard output files for each iterative step:

+-----------------------+------------------------------------------------------------------------+
| ``*.dist.new``        | distribution functions of the current step                             |
+-----------------------+------------------------------------------------------------------------+
| ``*.dpot.new``        | the final potential update, created by ``calc_update``                 |
+-----------------------+------------------------------------------------------------------------+
| ``*.dpot.<number>``   | for each postupdate script, the ``.dpot.new`` is saved and a new one   |
+-----------------------+------------------------------------------------------------------------+
|                       | is created                                                             |
+-----------------------+------------------------------------------------------------------------+
| ``*.pot.cur``         | the current potential used for the actual run                          |
+-----------------------+------------------------------------------------------------------------+
| ``*.pot.new``         | the new potential after the add step                                   |
+-----------------------+------------------------------------------------------------------------+
| ``*.pot.<number>``    | same as ``dpot.<number>`` but for ``post_add``                         |
+-----------------------+------------------------------------------------------------------------+

If a sub-step fails during the iteration, additional information can be
found in the log file. The name of the log file is specified in the
steering XMLfile.

Restarting and continuing
~~~~~~~~~~~~~~~~~~~~~~~~~

The interrupted or finished iterative process can be restarted either by
extending a finished run or by restarting the interrupted run. When the
script is called, it automatically checks for a file called ``done`` in
the current directory. If this file is found, the program assumes that
the run is finished. To extend the run, simply increase in the settings
file and remove the file called ``done``. After that, can be restarted,
which will automatically recognize existing steps and continue after the
last one.

If the iteration was interrupted, the script might not be able to
restart on its own. In this case, the easiest solution is to delete the
last step and start again. The script will then repeat the last step and
continue. However, this method is not always practical since sampling
and analysis might be time-consuming and the run might have only crashed
due to some inadequate post processing option. To avoid repeating the
entire run, the script creates a file with restart points and labels
already completed steps such as simulation, analysis, etc. The file name
is specified in the option . If specific actions should be redone, one
can simply remove the corresponding lines from this file. Note that a
file ``done`` is also created in each folder for those steps which have
been successfully finished.

Iterative Boltzmann Inversion
-----------------------------

Input preparation
~~~~~~~~~~~~~~~~~

This section describes the usage of IBI, implemented within the
scripting framework described in the previous section
[sec:iterative\_workflow]. It is suggested to get a basic understanding
of this framework before proceeding.

An outline of the workflow for performing IBIis given in
fig. [fig:flow\_ibi].

To specify Iterative Boltzmann Inversion as algorithm in the script, add
``ibi`` in the ``method`` section of the XMLsetting file as shown below.

::

      <cg>
        ...
        <inverse>
          <method>ibi</method>
        </inverse>
      </cg>

Inverse Monte Carlo
-------------------

In this section, additional options are described to run IMCcoarse
graining. The usage of IMCis similar to the one of IBIand understanding
the use of the scripting framework described in
chapter [sec:iterative\_workflow] is necessary.

**WARNING: multicomponent IMCis still experimental!**

General considerations
~~~~~~~~~~~~~~~~~~~~~~

In comparison to IBI, IMCneeds significantly more statistics to
calculate the potential update:raw-latex:`\cite{Ruehle:2009.a}`. It is
advisable to perform smoothing on the potential update. Smoothing can be
performed as described in sec. [ref:ibi:optimize]. In addition, IMCcan
lead to problems related to finite size: for methanol, an undersized
system proved to lead to a linear shift in the
potential:raw-latex:`\cite{Ruehle:2009.a}`. It is therefore always
necessary to check that the system size is sufficiently large and that
runlength csg smoothing iterations are well balanced.

Correlation groups
~~~~~~~~~~~~~~~~~~

Unlike IBI, IMCalso takes cross-correlations of interactions into
account in order to calculate the update. However, it might not always
be beneficial to evaluate cross-correlations of all pairs of
interactions. By specifying , votcaallows to define groups of
interactions, amongst which cross-correlations are taken into account,
where can be any name.

::

      <non-bonded>
        <name>CG-CG</name>
        <type1>CG</type1>
        <type2>CG</type2>
        ...
        <imc>
          <group>solvent</group>
       </imc>
      </non-bonded>
      <non-bonded>

Regularization
~~~~~~~~~~~~~~

To use the regularized version of IMC a :math:`\lambda` value :math:`>0`
has to be specified by setting . If set to :math:`0` (default value) the
unregularized version of IMC is applied.

::

     <non-bonded>
       <name>CG-CG</name>
       <type1>CG</type1>
       <type2>CG</type2>
        ...
       <inverse>
         <imc>
           <reg>300</reg>
         </imc>
       </inverse>
     </non-bonded>

Relative Entropy
----------------

In this section, additional options are described to run REcoarse
graining. The usage of REis similar to the one of IBIand IMCand
understanding the use of the scripting framework described in
chapter [sec:iterative\_workflow] is necessary.

Currently, REimplementation supports optimization of two-body non-bonded
pair interactions. Support for bonded and N-body interactions is
possible by further extension of REimplementation.

Potential function and parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In RE, CG potentials are modeled using analytical functional forms.
Therefore, for each CG interaction, an analytical functional must be
specified in the XMLsetting file as

::

      <non-bonded>
        <name>CG-CG</name>
        <type1>CG</type1>
        <type2>CG</type2>
        ...
        <re>
          <function>cbspl or lj126</function>
            <cbspl>
              <nknots>48</nknots>
            </cbspl>
        </re>
        ...
      </non-bonded>

Currently, standard Lennard-Jones 12-6 (lj126) and uniform cubic
B-splines-based piecewise polynomial (cbspl) functional forms are
supported. For lj126, the parameters to optimize are the usual
:math:`C_{12}` and :math:`C_{6}`. The cbspl form is defined as

.. math::

   \label{eq:cbspl}
   u_{\text{cbspl}}(r) = \left[\begin{array}{cccc}
       1 & t & t^2 & t^3 \end{array}\right]
   \frac{1}{6}
   \left[ \begin{array}{rrrr}
       1 & 4 & 1 & 0 \\
       -3 & 0 & 3 & 0 \\
       3 & -6 & 3 & 0 \\
       -1 & 3 & -3 & 1 \end{array}\right]
   \left[ \begin{array}{l}
       c_{k} \\
       c_{k+1} \\
       c_{k+2} \\
       c_{k+3} \end{array}\right] ,

 where :math:`\{c_0,c_1,c_2,...,c_m\}` are the spline knot values
tabulated for :math:`m` evenly spaced intervals of size
:math:`\Delta r = r_{\text{cut}}/(m-2)` along the separation distance
:math:`r_{i} = i\times\Delta r` with the cut-off :math:`r_{\text{cut}}`,
and :math:`t` is given by

.. math::

   \label{eq:cbspl_t}
   t = \frac{r-r_{k}}{\Delta r} ,

 where index :math:`k` is determined such that
:math:`r_{k}\leq r < r_{k+1}`. For cbspl, the knot values,
:math:`\{c_0,c_1,c_2,...,c_m\}`, are optimized. The number of knot
values to use must be specified in the XMLsetting file as shown in the
above snippet. :math:`u_{\text{cbspl}}(r)` exhibits remarkable
flexibility, and it can represent various complex functional
characteristics of pair potentials for sufficiently large number of
knots.

Update scaling parameter
~~~~~~~~~~~~~~~~~~~~~~~~

Depending on the quality of the initial guess and sensitivity of the CG
system to the CG parameters, scaling of the parameter update size may be
required to ensure the stability and convergence of the REminimization.
The scaling parameter, :math:`\chi\in(0...1)`, value can be specified in
the XMLsettings file.

Statistical averaging of parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Due to stochastic nature of the CG simulations, near convergence, the CG
potential paramters may fluctuate around the mean converged values.
Therefore, the optimal CG parameters can be estimated by averaging over
the last few iterations. To specify averaging, the ``average``, keyword
should be specified in the ``post_update`` options in the XMLsettings
file.

General considerations
~~~~~~~~~~~~~~~~~~~~~~

To ensure the stability of the relative entropy minimization, some
precautionary measures are taken. For the Newton-Raphson update to
converge towards a minimum, the Hessian, :math:`\mathbf{H}`, must be
positive definite at each step. With a good initial guess for the CG
parameters and by adjusting the value of the relaxation parameter,
:math:`\chi`, stability of the Newton-Raphson method can be ensured. One
approach to initialize the CG parameters can be to fit them to PMF
obtained by inverting the pair distributions of the CG sites obtained
from the reference AA ensemble. For the lj126 and cbspl forms, which are
linear in its parameters, the second derivative of
:math:`S_{\text{rel}}` is never negative, hence the minimization
converges to a single global minimum. However, due to locality property
of the cbspl form, i.e., update to :math:`c_i` affects only the value of
the potential near :math:`r_i`, and the poor sampling of the very small
separation distances in the high repulsive core, the rows of
:math:`\mathbf{H}` corresponding to the first few spline knots in the
repulsive core may become zero causing :math:`\mathbf{H}` to be a
singular matrix. To avoid this singularity issue, we specify a minimum
separation distance, :math:`r_{\text{min}}`, for each CG pair
interaction and remove the spline knots corresponding to the
:math:`r\le r_{\text{min}}` region from the Newton-Raphson update. Once
the remaining knot values are updated, the knot values in the poorly
sampled region, i.e., :math:`r\le r_{\text{min}}`, are linearly
extrapolated. The value of :math:`r_{\text{min}}` at each iteration is
estimated from the minimum distance at which the CG RDF from the CG-MD
simulation is nonzero. Also, to ensure that the CG pair potentials and
forces go smoothly to zero near :math:`r_{\text{cut}}`, 2 knot values
before and after :math:`r_{\text{cut}}`, i.e., total 4, are fixed to
zero.

Pressure correction
-------------------

The pressure of the coarse-grained system usually does not match the
pressure of the full atomistic system. This is because iterative
Boltzmann inversion only targets structural properties but not
thermodynamic properties. In order correct the pressure in such a way
that it matches the target pressure ()., different strategies have been
used based on small modifications of the potential. The correction can
be enable by adding pressure to the list of scripts. The type of
pressure correction is selected by setting .

Simple pressure correction
~~~~~~~~~~~~~~~~~~~~~~~~~~

In ref.:raw-latex:`\cite{Reith:2003}` a simple linear attractive
potential was added to the coarse-grained potential

.. math:: \Delta V(r)=A \left( 1-\frac{r}{r_{cutoff}} \right) \,,

 with prefactor :math:`A`

.. math:: A = -{\ensuremath{\operatorname{sgn}}}(\Delta P)0.1k_{B}T\min(1,|f\Delta P) \,,

 :math:`\Delta p=P_i-P_\text{target}`, and scaling factor :math:`f` and
:math:`P_\text{target}` can be specified in the settings file as and .

As an example for a block doing simple pressure correction, every third
interaction is

::

    <post_update>pressure</post_update>
    <post_update_options>
      <pressure>
        <type>simple</type>
        <do>0 0 1</do>
        <simple>
          <scale>0.0003</scale>
        </simple>
      </pressure
    </post_update_options>

Here, is the scaling factor :math:`f`. In order to get the correct
pressure it can become necessary to tune the scaling factor :math:`f`
during the iterative process.

Advanced pressure correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In :raw-latex:`\cite{Wang:2009}` a pressure correction based on the
virial expression of the pressure was introduced. The potential term
remains as in the simple form while a different sturcture of the
:math:`A` factor is used:

.. math:: A = \left[\frac{-2\pi\rho^{2}}{3r_{cut}}\int_{0}^{r_{cut}}r^{3}g_{i}(r)dr\right]A_{i}=\Delta P.

 This factor requires the particle density :math:` \rho ` as additional
input parameter, which is added as in the input file.

Kirkwood-Buff correction
------------------------

In order to reproduce the exact Kirkwood-Buff ingetrals (KBIs), an
correction term can be added into the coarse-grained
potential :raw-latex:`\cite{Ganguly:2012}`,

.. math:: \Delta U_{ij}^{(n)}(r) = \frac{k_{B}T}\;A\;(G_{ij}^{(n)} - G_{ij}^\text{ref})\left(1- \frac{r}{r_\text{ramp}}\right),

 where :math:`G_{ij}^{(ref)}` is the KBI calculated from the reference
all-atom simulation and :math:`G_{ij}^{(n)}` is the KBI after the
:math:`n^{th}` iteration.

The Kirkwood-Buff integrals are calculated from the radial distribution
functions as follows:

.. math::

   G_{ij} = 4\pi \int_0^\infty \left[ g_{ij}(r) - 1\right] r^2 dr~.
   \label{eq:kbi}

 For simulations of finite box size we calculate the running integral up
to distance :math:`R`

.. math:: G_{ij}(R) = 4\pi \int_0^R \left[ g_{ij}(r) - 1\right] r^2 dr~.

 The average of those running integrals in the interval, where
:math:`G_{ij}(R)` gets flat, gives a good estimate for :math:`G_{ij}`:

.. math:: G_{ij}\approx<G_{ij}(R)>|_{R=r_1}^{R=r_2}

 As an example for a block doing Kirkwood-Buff correction, every
iteraction without doing potential update

::

    <do_potential>0</do_potential>
    <post_update>kbibi</post_update>
    <post_update_options>
      <kbibi>
        <do>1</do>
        <start>1.0</start>
        <stop>1.4</stop>
        <factor>0.05</factor>
        <r_ramp>1.4</r_ramp>
      </kbibi>
    </post_update_options>

Here, is the scaling factor :math:`A`. is :math:`r_1` and is :math:`r_2`
used to calculate the average of :math:`G_{ij}(R)`.

Runtime optimization
--------------------

Most time per iteration is spent on running the coarse-grained system
and on calculating the statistics. To get a feeling on how much
statistics is needed, it is recommended to plot the distribution
functions and check whether they are sufficiently smooth. Bad statistics
lead to rough potential updates which might cause the iterative
refinement to fail. All runs should be long enough to produce
distributions/rdfs of reasonable quality.

Often, runtime can be improved by smoothing the potential updates. Our
experience has shown that it is better to smooth the potential update
instead of the rdf or potential itself. If the potential or rdf is
smoothed, sharp features like the first peak in SPC/Ewater might get
lost. Smoothing on the delta potential works quite well, since the sharp
features are already present from the initial guess. By applying
iterations of a simple triangular smoothing
(:math:` \Delta U_i = 0.25 \Delta U_{i-1} + 0.5\Delta U_i + 0.25\Delta U_{i+1} `),
a reasonable coarse-grained potential for SPC/Ewater could be produced
in less than 10 minutes. Smoothing is implemented as a post\_update
script and can be enabled by adding

::

      <post_update>smooth</post_update>
      <post_update_options>
        <smooth>
            <iterations>2</iterations>
        </smooth>
      </post_update_options>

to the inverse section of an interaction in the settings XMLfile.

Coordination Iterative Boltzmann Inversion
------------------------------------------

The method :math:`\mathcal{C}-`\ IBI (Coordination Iterative Boltzmann
Inversion) uses pair-wise cumulative coordination as a target function
within an iterative Boltzmann inversion. This method reproduces
solvation thermodynamics of binary and ternary mixtures
:raw-latex:`\cite{deOliveira:2016}`.

The estimation of coordination is given by:

.. math::

   \label{eq:coord}
   \mathcal{C}_{ij}(r) = 4\pi \int_{0}^{r} {\rm g}_{ij}(r')r'^{2}dr'

 with the indices :math:`i` and :math:`j` standing for every set of
pairs, uses a volume integral of :math:`{\rm g}(r)`.

The Kirkwood and Buff theory (KB) :raw-latex:`\cite{Kirkwood:1951}`
connects the pair-wise coordinations with particule fluctuations and,
thus, with the solution thermodynamics
:raw-latex:`\cite{Mukherji:2013,Naim:2006}`. This theory make use of the
Kirkwood-Buff integrals (KBI) :math:`{\rm G}_{ij}` defined as,

.. math::

   \label{eq:Gij}
   {\rm G}_{ij} = 4 \pi \int_{0}^{\infty} \left [ {\rm g}_{ij}(r) - 1 \right ] r^{2} dr.

 For big system sizes the :math:`{\rm G}_{ij}` can be approximated:

.. math::

   \label{eq:Gij_app}
   {\rm G}_{ij} = \mathcal{C}_{ij}(r) - \frac{4}{3} \pi r^{3},

 were the second therm is a volume correction to
:math:`\mathcal{C}_{ij}(r)`.

Thus the initial guess for the potential of the CG model is obtained
from the all atom simulations,

.. math::

   \label{eq:pot_ibi}
   {\rm V}_{0}(r) = -k_{B}T {\rm ln} \left [ {\rm g}_{ij}(r) \right ],

 however, the iterative protocol is modified to target
:math:`\mathcal{C}_{ij}(r)` given by,

.. math::

   \label{eq:pot_cibi}
   {\rm V}_{n}^{\mathcal{C}-{\rm IBI}}(r) = {\rm V}_{n-1}^{\mathcal{C}-{\rm IBI}}(r)
   + k_{B}T {\rm ln} \left [ \frac{\mathcal{C}_{ij}^{n-1}(r)}{\mathcal{C}_{ij}^{target}(r)} \right ].

To perform the :math:`\mathcal{C}-`\ IBI is necessary include some lines
inside of the .xml file:

::

     <cg>
      <non-bonded>
       <name>A-A</name>
       ...
       <inverse>
        <post_update>cibi</post_update>
        <post_update_options>
          <cibi>
            <do>1</do>
          </cibi>
        </post_update_options>
        ...
     </cg>

DL\_POLYinterface
=================

**WARNING: The DL\_POLYinterface is still experimental (in development)
but it does support the Iterative Boltzmann Inversion and Inverse Monte
Carlo schemes. The Force Matching might work as well, although it has
not been tested thoroughly.**

General remarks on using  votcawith DL\_POLY
-------------------------------------------

The DL\_POLYinterface fully supports coarse-grain mapping of a full-atom
system previuosly simulated with any version of DL\_POLY, including
DL\_POLY-Classic. However, the full optimization of the effective
potentials with the aid of iterative methods will only become possible
when the new release of DL\_POLY-4(4.06) is made public; the reason
being the incapability of earlier DL\_POLYversions of using
user-specified tabulated force-fields for intramolecular, aka “bonded”,
interactions: bonds, angles, dihedral angles (torsions). Below the
coarse-graining and CG force-field optimization with the aid of the
latest DL\_POLY-4version (4.06+) are outlined.

Running votcawith DL\_POLY-4as MD simulation engine is very similar to
doing so with GROMACS. The three types of required input files in the
case of DL\_POLYare: CONTROL – containing the simulation directives and
parameters (instead of ``.mdp`` file for GROMACS), FIELD – the topology
and force-field specifications (instead of ``.top`` and ``.tpr`` files),
and CONFIG (instead of ``.gro`` file) – the initial configuration file,
containing the MD cell matrix and particle coordinates (it can also
include initial velocities and/or forces); for details see
DL\_POLY-4manual. Most of the votcatools and scripts described above in
the case of using GROMACSwill work in the same manner, with the
following conventional substitutions for the (default) file names used
in options for votcascripts, as necessary:

::

    .dlpf = the topology read from FIELD or written to FIELD_CGV
    .dlpc = the configuration read from CONFIG or written to CONFIG_CGV
    .dlph = the trajectory read from HISTORY or written to HISTORY_CGV

It is also possible to specify file names different from the standard
DL\_POLYconvention, in which case the user has to use the corresponding
dot-preceded extension(s); for example: FA-FIELD.dlpf instead of FIELD
or CG-HISTORY.dlph instead of HISTORY\_CGV (see section
[sec:ref\_programs], as well as the man pages or output of
votcacommands, with option ``—help``).

 votcafollows the DL\_POLYconventions for file names and formats. Thus,
``csg_dlptopol`` and ``csg_map`` produce the CG topology (FIELD\_CGV by
default), configuration (CONFIG\_CGV), and/or trajectory (HISTORY\_CGV)
files fully compatible with and usable by DL\_POLY. **Note that the
ability of these tools to read and write a plethora of different file
formats provides means to convert input and output files between the
simulation packages supported by votca, e.g. GROMACS– DL\_POLYor vice
versa. The user is, however, strongly advised to check the resulting
files for consistency before using them).**

Similarly, the distribution analysis and potential/force generation
utilities, such as ``csg_stat`` and votcascripts, will read and write
DL\_POLY-formatted files; in particular, the tabulated force-field files
containing the potential and force/virial data: TABLE – for short-range
(VdW) “non-bonded” interactions, TABBND, TABANG and TABDIH – for
“bonded” interations: bonds, bending angles and dihedrals,
correspondingly (for the format details see DL\_POLY-4manual). Note,
however, that the latter three files can only be used by
DL\_POLY-4(4.06+).

The user is advised to search for “dlpoly” through the
``csg_defaults.xml``, ``csg_table`` files and in scripts located in
``share/votca/scripts/inverse/`` in order to find out about the xml-tags
and options specific for DL\_POLY; see also sections [sec:ref\_options]
and [sec:csg\_table].

Advanced topics
===============

Customization
-------------

Each sub-step of an iteration and all direct calls can be adjusted to
the user needs. The internal part of the iterative framework is
organized as follows: all scripts are called using two keywords

::

      csg_call key1 key2

For example, ``csg_call update imc`` calls the ``update`` script for the
inverse Monte Carlo procedure. The corresponding keywords are listed in
sec. [sec:csg\_table] or can be output directly by calling

::

      csg_call --list

It is advised not to change already implemented scripts. To customize a
script or add a new one, copy the script to your own directory (set by )
and redirect its call by creating your own ``csg_table`` file in this
directory which looks like this

::

      key1 key2 script1 options
      key3 key4 script2

If the local keys are already in use, the existing call will be
overloaded.

As an example, we will illustrate how to overload the script which calls
the sampling package. The script runs from the GROMACSpackage only on
one cpu. Our task will be to change the script so that GROMACSuses 8
cpus, which is basically the same as adding mpirun options in .

First we find out which script calls :

::

      csg_call --list | grep gromacs

The output should look as follows

::

      init gromacs initalize_gromacs.sh
      prepare gromacs prepare_gromacs.sh
      run gromacs run_gromacs.sh
      pressure gromacs calc_pressure_gromacs.sh
      rdf gromacs calc_rdf_gromacs.sh
      imc_stat gromacs imc_stat_generic.sh
      convert_potential gromacs potential_to_gromacs.sh

the third line indicates the script we need. If the output of is not
clear, one can try to find the right script in sec. [sec:csg\_table].
Alternatively, check the folder

::

      <csg-installation>/share/scripts/inverse

for all available scripts.

Analyzing the output of

::

      csg_call --cat run gromacs

we can conclude that this is indeed the script we need as the content
(in shorted form is):

::

      critical mdrun

Now we can create our own ``SCRIPTDIR``, add a new script there, make it
executable and overload the call of the script:

::

      mkdir -p SCRIPTDIR
      cp `csg_call --quiet --show run gromacs` SCRIPTDIR/my_run_gromacs.sh
      chmod 755 SCRIPTDIR/my_run_gromacs.sh
      echo "run gromacs my_run_gromacs.sh" >> SCRIPTDIR/csg_table

Please note that ``my_run_gromacs.sh`` is the name of the script and
``SCRIPTDIR`` is the custom script directory, which can be a global or a
local path. Now we change the last line of ``my_run_gromacs.sh`` to:

::

      critical mpirun -np 8 mdrun

This completes the customization. Do not forget to add ``SCRIPTDIR`` to
in the setting XMLfile (see sec. [sec:ref\_options]).

You can check the new script by running:

::

      csg_call --scriptdir SCRIPTDIR --list
      csg_call --scriptdir SCRIPTDIR --run run gromacs

Finally, do not forget to remove the license infomation and change the
version number of the script.

Used external packages
----------------------

GroMaCS
~~~~~~~

Get it from

-  
-  

ESPResSo
~~~~~~~~

Get it from

DL\_POLY
~~~~~~~~

Get it from

Gnuplot
~~~~~~~

Get it from

LAMMPS
~~~~~~

Get it from

Reference
=========

Programs
--------

Mapping file
------------

The root node always has to be cg\_molecule. It can contain the
following keywords:

Topology file
-------------

The XMLtopology file

Settings file
-------------

All options for the iterative script are stored in an xml file.
[sec:ref\_options]

Scripts
-------

Scripts are used by and . The script table commonly used (compare
``csg_call –list``):

Script calls can be overwritten by adding a line with the 3rd column
changed to ``csg_table`` in directory.

.. [1]
   In what follows we adopt notations of
   ref. :raw-latex:`\cite{Noid:2008.1}`.

.. [2]
   Checking the linear correlation coefficient does not guarantee
   statistical independence of variables, for example
   :math:`c(x, x^2)=0` if :math:`x` has a symmetric probability density
   :math:`P(x) = P(-x)`. This case is often encountered in systems used
   for coarse-graining.

.. |image| image:: fig/propane
   :width: 4.00000cm
.. |image| image:: functionality/fig/iterativemethods.eps
   :width: 3.00000cm
.. |image| image:: usage/fig/propane.eps
   :width: 6.00000cm
.. |image| image:: usage/fig/flow_boltzmann.eps
   :width: 7.00000cm
