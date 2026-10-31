Introduction
############

VOTCA-XTP calculates electronic excitations of molecules and molecular
materials and the transport of charges and excitons between them:

-  **Electronic structure** of single molecules and clusters: density-functional
   theory (the internal DFT engine or ORCA), quasiparticle energies in the GW
   approximation and neutral excitations from the Bethe-Salpeter equation
   (GW-BSE), see :doc:`elstruc_theory`, :doc:`implementation_notes` and
   :doc:`gwbse_numerics`.
-  **Excitations in a molecular environment**: quantum-classical (QM/MM) embedding
   with polarizable classical regions (Thole or charge-response models,
   cut-off or Ewald summation) and embedded GW-BSE, see
   :doc:`ewald_embedding`, :doc:`polarizable_environment` and
   :doc:`embedded_gwbse`.
-  **Charge and exciton transport** in morphologies from molecular dynamics:
   site energies and electronic couplings for all neighbouring molecules,
   Marcus rates, kinetic Monte Carlo, time-of-flight and graph random walks,
   and stochastic network models for large systems, see
   :doc:`transfer_theory`, :doc:`transport_theory` and
   :doc:`transport_simulations`.

The programs are

-  ``xtp_tools``: calculations that do not need a morphology (e.g. ``dftgwbse``
   for a single molecule),
-  ``xtp_map``: maps an MD trajectory to a state file (``state.hdf5``) with
   the molecules ("segments") of the morphology,
-  ``xtp_run``: calculators that act on the whole state file (neighbour list,
   analysis, KMC and other transport calculations),
-  ``xtp_parallel``: job-based calculators that run many independent
   quantum-chemical or QM/MM calculations on the state file.

Their options are listed in :doc:`workflows_and_calculators`; ``-d <name>``
prints the options of a calculator and ``--cite <name>`` the papers to cite
for it. A Python (ASE) interface is described in :doc:`pyxtp`. The tutorials
show complete workflows.
