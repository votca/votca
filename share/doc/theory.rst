Theoretical background
======================

.. _theory_mapping:

Mapping
-------

The mapping is an operator that establishes a link between the atomistic
and coarse-grained representations of the system. An atomistic system is
described by specifying the values of the Cartesian coordinates and
momenta

.. math::

   \begin{aligned}
   \mathbf r^n &= \{\mathbf r_1,\dots,\mathbf r_n\}, \\
   \mathbf p^n &= \{\mathbf p_1,\dots,\mathbf p_n\}.
   \end{aligned}

of the :math:`n` atoms in the system. [#f1]_ On a coarse-grained level,
the coordinates and momenta are specified by the positions and momenta
of CG sites

.. math::

   \begin{aligned}
   \mathbf R^N = \{\mathbf R_1,\dots,\mathbf R_N\}, \\
   \mathbf P^N = \{\mathbf P_1,\dots,\mathbf P_N\}.
   \end{aligned}

Note that capitalized symbols are used for the CG sites while lower
case letters are used for the atomistic system.

The mapping operator :math:`{\mathbf c}_I` is defined by a matrix for each
bead :math:`I` and links the two descriptions

.. _theory_eq_mapping_scheme:

.. math::

   \begin{aligned}
    {\mathbf R}_I &= \sum_{i=1}^{n}c_{Ii}\mathbf r_i, \\
    {\mathbf P}_I &=
       M_I \dot{{\mathbf R}}_I =
       M_I \sum_{i=1}^{n}c_{Ii} \dot{{\mathbf r}}_i =
       M_I \sum_{i=1}^{n} \frac{ c_{Ii}} {m_i} {\mathbf p}_i .
   \end{aligned}

for all :math:`I = 1,\dots,N`.

If an atomistic system is translated by a constant vector, the
corresponding coarse-grained system is also translated by the same
vector. This implies that, for all :math:`I`,

.. math:: \sum_{i=1}^{n}c_{Ii}=1.

In some cases it is useful to define the CG mapping in such a way that
certain atoms belong to several CG beads at the same
time [Fritz:2009]_. Following
ref. [Noid:2008.1]_, we define two sets of atoms for
each of the :math:`N` CG beads. For each site :math:`I`, a set of
*involved* atoms is defined as

.. math:: {\cal I}_I=\{i|c_{Ii}\ne0\}.

An atom :math:`i` in the atomistic model is involved in a CG site, *I*,
if and only if this atom provides a nonzero contribution to the sum in
:ref:`the equation above<theory_eq_mapping_scheme>`.

A set of *specific* atoms is defined as

.. math:: {\cal S}_I=\{i|c_{Ii}\ne0 \text{ and } c_{Ji}=0 \text{ for all } J \ne I\}.

In other words, atom :math:`i` is specific to site :math:`I` if and
only if this atom is involved in site :math:`I` and is not involved in
the definition of any other site.

The CG model will generate an equilibrium distribution of momenta that
is consistent with an underlying atomistic model if all the atoms are
*specific* and if the mass of the :math:`I^\text{th}` CG site is given
by [Noid:2008.1]_

.. math::

   M_I= \left( \sum_{i \in {\cal I}_I}\frac{c_{Ii}^2}{m_i} \right)^{-1}.

If all atoms are specific and the center of mass of a bead is used for
mapping, then :math:`c_{Ii} = \frac{m_i}{M_I}`, and the
condition is automatically satisfied.

.. rubric:: Footnotes
.. [#] In what follows we adopt notations of ref. [Noid:2008.1]_.



Boltzmann inversion
-------------------

Boltzmann inversion is mostly used for *bonded* potentials, such as
bonds, angles, and torsions [Tschoep:1998]_. Boltzmann
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

.. _theory_eq_boltzmann_norm:

.. math::

   \begin{aligned}
   P_r(r) = \frac{H_r(r)}{4\pi r^2}~,\;
   P_\theta(\theta) = \frac{H_\theta(\theta)}{\sin \theta}~,\;
   P_\varphi(\varphi) = H_\varphi (\varphi)~,
   \end{aligned}

where :math:`r` is the bond length :math:`r`, :math:`\theta` is the
bond angle, and :math:`\varphi` is the torsion angle. The bonded
coarse-grained potential can then be written as a sum of distribution
functions

.. _theory_eq_boltzmann_pmf:

.. math::

   \begin{aligned}
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
probability distribution and the potential :ref:`as in the above equation<theory_eq_boltzmann_pmf>`.
Hence, one has to carefully check whether this assumption holds in
practice. This can be done by performing coarse-grained simulations and
comparing cross-correlations for all pairs of degrees of freedom in
atomistic and coarse-grained resolution, e. g. using a two-dimensional
histogram, analogous to a Ramachandran plot.  [#f2]_

.. rubric:: Footnotes
.. [#] Checking the linear correlation coefficient does not guarantee
   statistical independence of variables, for example
   :math:`c(x, x^2)=0` if :math:`x` has a symmetric probability density
   :math:`P(x) = P(-x)`. This case is often encountered in systems used
   for coarse-graining.

Separation of bonded and non-bonded degrees of freedom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When coarse-graining polymeric systems, it is convenient to treat bonded
and non-bonded interactions separately [Tschoep:1998]_.
In this case, sampling of the atomistic system shall be performed on a
special system where non-bonded interactions are artificially removed,
so that the non-bonded interactions in the reference system do not
contribute to the bonded interactions of the coarse-grained model.

This can be done by employing exclusion lists using with the option
``—excl``. This is described in detail in :ref:`methods_exclusions`.

.. figure:: fig/excl.png
   :align: center

   Example of excluded interactions.

.. _theory_iterative_methods:

Iterative methods
-----------------

.. _theory_fig_iterative_simple:

.. figure:: fig/iteration-scheme-simple.png
   :align: center

   Block-scheme of an iterative method.

Iterative workflow control is essential for the IBI and IMC methods. The
general idea of iterative workflow is sketched in
:ref:`the block-scheme above<theory_fig_iterative_simple>`. A run starts with
an initial guess during the global initialization phase. This guess is used for
the first sampling step, followed by an update of the potential. The update
itself often requires additional postprocessing such as smoothing,
interpolation, extrapolation or fitting. Different methods are available to
update the potential, for instance Iterative Boltzmann Inversion (see
:ref:`theory_iterative_boltzmann_inversion`) or Inverse Monte Carlo (see
:ref:`theory_inverse_monte_carlo`). The whole procedure is then iterated until
a convergence criterion is satisfied.

.. _theory_iterative_boltzmann_inversion:

Iterative Boltzmann Inversion
-----------------------------

Iterative Boltzmann inversion (IBI) is a natural extension of the
Boltzmann inversion method. Since the goal of the coarse-grained model
is to reproduce the distribution functions of the reference system as
accurately as possible, one can also iteratively refine the
coarse-grained potentials using some numerical scheme.

In IBIthe potential update :math:`\Delta U` is given
by [Reith:2003]_

.. math::

    \begin{aligned}
    U^{(n+1)} &= U^{(n)} + \lambda \Delta U^{(n)}~, \\
    \Delta U^{(n)} &=  k_\text{B} T \ln  \frac{P^{(n)}}{P_{\rm ref}}
    =  U_\text{PMF}^\text{ref} - U_\text{PMF}^{(n)}~.
    \end{aligned}
    \label{eq:iter_boltzmann}

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

.. _theory_inverse_monte_carlo:

Inverse Monte Carlo
-------------------

Inverse Monte Carlo (IMC) is an iterative scheme which additionally
includes cross correlations of distributions. A detailed derivation of
the IMCmethod can be found in ref. [Lyubartsev:1995]_.

The potential update :math:`\Delta U` of the IMCmethod is calculated by
solving a set of linear equations

.. math::

   \begin{aligned}
       \left<S_{\alpha}\right> - S_{\alpha}^{\text{ref}}= A_{\alpha \gamma} \Delta U_{\gamma}~,
     \label{eq:imc}\end{aligned}

where

.. _theory_eq_covariance:

.. math::

   \begin{aligned}
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
reformulate :ref:`the above equation<theory_eq_covariance>` by addition of a penalty term. In
this case the potential update is computed as
follows:[Murtola:2007]_

.. math::

   \Delta U_\gamma = \arg \min \| A_{\alpha \gamma} \Delta U_\gamma - \left(\left<S_{\alpha}\right> - S_{\alpha}^{\text{ref}}\right) \|^2 + \lambda \| R \Delta U_{\gamma} \|^{2}

This equation is known as Tikhonov regularization, where
:math:`R` is the regularization operator, which here is the identity
matrix and :math:`\lambda >0 ` is the regularization parameter. The
optimal choice for :math:`\lambda` can only be determined if the exact
solution of the inverse problem is known, which in practice is not the
case. To get a good initial guess on the magnitude of the regularization
parameter a singular value decomposition of the matrix
:math:`A_{\alpha \gamma}` might help. A good :math:`\lambda` parameter
should dominate the smallest singular values (squared) but is itself
small compared to the larger ones.[Rosenberger:2016]_

Iterative Integral Equation methods
-----------------------------------

Iterative Integral Equation (IIE) methods are best compared to the IMC method.
The main difference is that the Jacobian is not sampled from particle number
fluctuations, but approximately infered from the RDFs. The connection of
potential and RDF is obtained from inverting the Ornstein-Zernicke (OZ)
equation and a closing relation, e.g. the hypernetted-chain (HNC) equation

.. math::

   \begin{aligned}
   h &= c + \rho h * c \\
   u &= k_B T \left(h - \ln(g) - c\right)
   \end{aligned}

where :math:`h` is :math:`g - 1`,
:math:`c` is the direct correlation function,
:math:`\rho` is the density,
:math:`u` is the pair potential,
and :math:`*` denotes a 3D convolution.
For the case of bonds, the reference interaction site model (RISM) form of the
OZ equation is used. For multiple bead types, the OZ equation becomes a matrix
equation but this case is currently not implemented.

The Gauss-Newton formalism allows the incorporation of one or multiple
constraints into the potential update. Those have to be expressible in terms
of the RDF and the potential or force. Currently, only a pressure constraint is
implemented, which is defined by

.. math::

   \begin{equation}
   \vec{l^T} \vec{f_k} = p_\text{tgt} - p_k .
   \end{equation}

Here :math:`p_\text{tgt}` and :math:`p_k` are the target and current pressure,
respectively, and :math:`\vec{f_k}` is the current pair force.
Element :math:`\alpha` of vector :math:`\vec{l}` is defined as

.. math::

    \begin{equation}
    l_\alpha = \frac{2}{3} \pi \rho^2 \frac{g_{\text{tgt},\alpha} + g_{\text{tgt},\alpha+1}}{2} \frac{r_{\alpha+1}^4 - r_\alpha^4}{4}
    \end{equation}

where :math:`\rho` is the particle density and :math:`r` is the radius.

The exact formulas and their derivation can be found in ref. [Delbary:2020]_
and [Bernhardt:2021]_. 

Force Matching
--------------

Force matching (FM) is another approach to evaluate corse-grained
potentials [Ercolessi:1994,Izvekov:2005,Noid:2007]_. In
contrast to the structure-based approaches, its aim is not to reproduce
various distribution functions, but instead to match the multibody
potential of mean force as close as possible with a given set of
coarse-grained interactions.

The method works as follows. We first assume that the coarse-grained
force-field (and hence the forces) depends on :math:`M` parameters
:math:`g_1,...,g_M`. These parameters can be prefactors of analytical
functions, tabulated values of the interaction potentials, or
coefficients of splines used to describe these potentials.

In order to determine these parameters, the reference forces on
coarse-grained beads are calculated by summing up the forces on the
atoms

.. math::

   {{{{\mathbf F}}}}_I^\text{ref} = \sum_{j \in {\cal S_I}} \frac{d_{Ii}}{c_{Ii}} {{{{\mathbf f}}}}_j({{{{\mathbf r}}}^n}),
     \label{eq:force_mapping}

where the sum is over all atoms of the CG site *I* (see
:ref:`theory_mapping`). The :math:`d_{Ij}` coefficients can, in
principle, be chosen arbitrarily, provided that the condition
:math:` \sum_{i=1}^{n}d_{Ii}=1` is
satisfied [Noid:2008.1]_. If mapping coefficients for
the forces are not provided, it is assumed that :math:`d_{Ij} = c_{Ij}`
(see also :ref:`input_files`).

By calculating the reference forces for :math:`L` snapshots we can write
down :math:`N \times L` equations

.. _theory_eq_fmatch1:

.. math::

   {{{{\mathbf F}}}}_{Il}^\text{cg}(g_1, \dots ,g_M)={{{\mathbf F}}}_{il}^\text{ref},\;
     I=1,\dots,N,\; l=1,\dots,L~.

Here :math:`{{{{\mathbf F}}}}_{Il}^\text{ref}` is the force on
the bead :math:`I` and :math:`{{{{\mathbf F}}}}_{Il}^\text{cg} `
is the coarse-grained representation of this force. The index :math:`l`
enumerates snapshots picked for coarse-graining. By running the
simulations long enough one can always ensure that
:math:`M < N times L`. In this case the set of equations
is overdetermined and can be solved in a least-squares manner.

:math:`{\mathbf F}_{il}^\text{cg}` is, in principle, a non-linear function
of its parameters :math:`\{g_i\}`. Therefore, it is useful to represent
the coarse-grained force-field in such a way that
:ref:`equations <theory_eq_fmatch1>` become linear functions of :math:`{g_i}`.
This can be done using splines to describe the functional form of the
forces [Izvekov:2005]_. Implementation details are
discussed in ref. [Ruehle:2009.a]_.

Note that an adequate sampling of the system requires a large number of
snapshots :math:`L`. Hence, the applicability of the method is often
constrained by the amount of memory available. To remedy the situation,
one can split the trajectory into blocks, find the coarse-grained
potential for each block and then perform averaging over all blocks.

Relative Entropy
----------------

Relative entropy is a method which quantifies the extent of the
configurational phase-space overlap between two molecular
ensembles [Wu2005]_. It can be used as a measure of the
discrepancies between various properties of the CG system’s and the
target all-atom (AA) ensemble. It has been shown by Shell
S. [Shell2008]_ that one can minimize the relative
entropy metric between the model CG system and the target AA system to
optimize CG potential parameters such that the CG ensemble would mimic
the target AA ensemble.

Relative entropy, :math:`S_{\text{rel}}`, is defined as
[Shell2008]_

.. _theory_eq_srel:

.. math::

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
probabilities into :ref:`the definition for the relative entropy<theory_eq_srel>`,
it simplifies to

.. _theory_eq_srelcan:

.. math::

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
ref. [Chaimovich2011]_. In this strategy, the CG
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
:math:`boldsymbollambda`, which can be computed from `the above equation<theory_eq_srelcan>`
as

.. math::

   \nabla_{\lambda}S_{\text{rel}} = \beta \left\langle \frac{\partial
     U_{\text{CG}}}{\partial\lambda}\right\rangle_{\text{AA}} - \beta\left\langle
   \frac{\partial U_{\text{CG}}}{\partial\lambda}\right\rangle_{\text{CG}} ,

and :math:`\mathbf{H}` is the Hessian matrix of :math:`S_{\text{rel}}`
given by

.. math::

   \begin{aligned}
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
from those two equations, we need average CG energy
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
are made in ref. [rudzinski_coarse-graining_2011]_
and [Chaimovich2011]_. Chaimovich and
Shell [Chaimovich2011]_ have shown that for certain CG
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
al [lyubartsev2010systematic]_ have shows how to use
IMC with an analytical function form, too. BI, IBI, and IMC methods are
based on pair correlations and hence, they are only useful to optimize
2-body CG potentials, whereas, relative entropy uses more generic metric
which offers more flexibility in modeling CG interactions and not only
2-body, but also 3-body (for example see
ref. [lu_coarse-graining_2014]_) and N-body CG
potentials can be optimized. In addition to the CG potential
optimization, the relative entropy metric can also be used to optimize
an AA to CG mapping operator.


