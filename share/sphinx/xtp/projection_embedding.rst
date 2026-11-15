Projection-Based Embedding (DFT-in-DFT)
#######################################

Many systems of interest are too large for a GW-BSE calculation of the
whole, while the excitation of interest is confined to a part of it: a
chromophore with long side chains, a molecule at a donor-acceptor
interface, a dye in its solvation shell. Projection-based embedding
[Manby:2012]_ splits such a system into an *active* region A, which keeps
its full electronic structure, and an *inactive* region B, which is
frozen at the DFT level. At the DFT level the split is exact: with the
same functional in both regions, the embedded calculation reproduces the
total energy and density of the full system. GW-BSE then works with the
orbitals of the active region only, which reduces the size of the
three-centre integrals and of the BSE Hamiltonian [Sundaram:2024]_.

Theory
******

The starting point is a Kohn-Sham calculation of the full system with
density matrix :math:`\gamma`. Its occupied orbitals are localized with
the Pipek-Mezey scheme [Pipek:1989]_. A localized orbital belongs to the
active region if its Mulliken population on at least one of the active
atoms exceeds a threshold (0.4 by default); all others belong to B. This
splits the density,

.. math::

   \gamma = \gamma_A^0 + \gamma_B .

The active region is then optimized in the field of the frozen inactive
density. With :math:`g[\gamma]` the Coulomb and exchange-correlation
potential of a density and the embedding potential

.. math::

   v_\text{emb} = g[\gamma_A^0 + \gamma_B] - g[\gamma_A^0] ,

the active electrons see the Fock matrix

.. math::

   F_A = h + g[\gamma_A] + v_\text{emb} + \mu P_B ,
   \qquad P_B = S \gamma_B S ,

and the total energy is

.. math::

   E[\gamma_A] = E_\text{DFT}[\gamma_A] + E_\text{DFT}[\gamma]
   - E_\text{DFT}[\gamma_A^0]
   + \mathrm{tr}\left[(\gamma_A - \gamma_A^0)\, v_\text{emb}\right]
   + \mu\, \mathrm{tr}\left[\gamma_A P_B\right] .

The projector :math:`P_B` with a large level shift :math:`\mu` keeps the
active orbitals orthogonal to the inactive ones: in the embedded
calculation the inactive orbitals appear as levels near :math:`\mu`, far
above all physical levels, and the active electrons cannot fall into
them. At self-consistency :math:`\gamma_A = \gamma_A^0`, and the energy
equals the full one; the log reports the difference (typically
:math:`10^{-7}` Ha).

The embedded calculation can also run in a basis restricted to the
active atoms and the *border* atoms of the inactive region, i.e. those
with a basis function whose Mulliken population in :math:`\gamma_A^0`
exceeds ``truncation_threshold`` [Barnes:2013]_, [Bennie:2015]_. This
makes the embedded calculation itself cheaper, but the result then
depends on the threshold; [Sundaram:2024]_ found it necessary to tighten
it from :math:`10^{-4}` to :math:`10^{-5}` for a charge-transfer complex.

Running an embedded calculation
*******************************

In ``xtp_tools -e dftgwbse`` the embedding needs the tasks ``localize``
and ``dft_in_dft`` after the DFT calculation of the full system, and the
active atoms in ``dftpackage.xtpdft.dft_in_dft``:

.. code-block:: xml

   <dftgwbse>
     <job_name>system</job_name>
     <dftpackage>
       <basisset>def2-svp</basisset>
       <auxbasisset>aux-def2-svp</auxbasisset>
       <functional>XC_HYB_GGA_XC_PBEH</functional>
       <xtpdft>
         <dft_in_dft>
           <activeatoms>0:15 105 106</activeatoms>
           <threshold>0.4</threshold>
           <levelshift>100000</levelshift>
           <truncate_basis>false</truncate_basis>
         </dft_in_dft>
       </xtpdft>
     </dftpackage>
     <tasks>input,dft,parse,localize,dft_in_dft,gwbse</tasks>
     <localize>
       <method>jacobi-sweeps</method>
     </localize>
     <gwbse>
       ...
     </gwbse>
   </dftgwbse>

The atom indices refer to the order of the atoms in the input geometry
and use the range syntax ``a:b`` (inclusive) with single indices,
separated by spaces or commas.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - option (``dft_in_dft``)
     - meaning
   * - ``activeatoms``
     - atoms of the active region
   * - ``threshold``
     - Mulliken population on an active atom above which a localized
       orbital is active (default 0.4)
   * - ``levelshift``
     - :math:`\mu` in Hartree (default :math:`10^4`)
   * - ``truncate_basis``, ``truncation_threshold``
     - restrict the basis to active and border atoms (default off,
       :math:`10^{-4}`)
   * - ``grow``, ``grow_window``, ``grow_threshold``
     - grow the active region from the frontier orbitals (see below;
       default off, 0.1 eV, 0.0025)

The localization (``localize.method``) is done with a unitary optimizer
by default; ``jacobi-sweeps`` can be more robust for large
molecules.

In a ``qmmm`` job the same ``dft_in_dft`` block goes into the
``dftpackage`` of the QM region, together with a ``localize`` block;
the active atoms are then numbered in the order of the atoms of the QM
region. The classical regions see the full density of the QM region.

Choosing the active region
**************************

The active region is a modelling choice, and the results of an embedded
GW-BSE calculation depend on it more than on anything else. Two things
should hold:

* the orbitals that take part in the excitation lie in it: an exciton
  whose hole extends across the cut cannot be described by the active
  orbitals;
* the cut goes through bonds that are electronically unimportant, such
  as the single bond to a saturated side chain, and not through a
  conjugated system.

To help with this, the embedding step reports, before it starts:

* whether the active atoms of each molecule form one connected piece of
  the bond graph (a region that splits a molecule stops the run);
* the coverage of the frontier orbitals of the full calculation: the
  fraction of each occupied orbital within ``grow_window`` of the HOMO
  (Mulliken population) that lies on the active region, and the inactive
  atoms on which they have weight;
* inactive hydrogens bonded to active atoms and inactive heavy atoms
  bonded to two or more active ones.

With ``grow`` the region is extended automatically: every atom of the
molecules in the active region whose population in one of these frontier
orbitals exceeds ``grow_threshold`` is added, together with the shortest
chain of bonds that connects it to the region, then heavy atoms bonded
to two or more active atoms, and the hydrogens of all active atoms.
Without a truncated basis only the occupied frontier orbitals matter,
since the embedded calculation keeps all virtual orbitals; with
``truncate_basis`` the virtual frontier orbitals are included as well.

GW-BSE in the active region
***************************

GW-BSE uses the occupied orbitals of the active region and, without a
truncated basis, all virtual orbitals of the full basis; the inactive
levels shifted to :math:`\mu` are left out. The QP and BSE windows refer
to these orbitals. A frozen core (``ignore_corelevels``) counts only the
core levels of the active atoms.

Two effects of the inactive region are missing in this plain form: it
does not screen the interaction (the QP gap comes out too large and the
electron-hole attraction too strong; the two errors partly cancel in the
excitation energies), and its occupied orbitals are missing from the
Green's function. Both can be added with ``gwbse.inactive_screening``
and ``gwbse.inactive_screening_complete_g``, described in
:ref:`pbe-inactive-screening`. The active region can be embedded in a
polarizable environment at the same time (:doc:`embedded_gwbse`).

The embedded GW-BSE is implemented for closed-shell active regions with
an even number of electrons.
