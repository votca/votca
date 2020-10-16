Coordination Iterative Boltzmann Inversion
==========================================

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
