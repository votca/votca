**********
XTP Manual
**********

**Disclaimer:**
This manual is not complete. The best way to start using the software is
to look at the provided tutorials. The reference section is generated
automatically from the source code of the master version, so please make sure that your
software and manual versions match.


.. toctree::
   :numbered:

   xtp_introduction
   elstruc_theory
   implementation_notes
   gwbse_numerics
   transfer_theory
   transport_theory
   transport_simulations
   ewald_embedding
   polarizable_environment
   embedded_gwbse
   pyxtp
   workflows_and_calculators
   reference
   bibliography


Citations
---------

Development of this software depends on academic research grants. If
you are using VOTCA-XTP, please cite [Baumeier:2024]_ and the papers
of the methods you used.

Every run of ``xtp_run``, ``xtp_tools`` and ``xtp_parallel`` ends with
the list of references for the methods it actually used, in the order
of first use, e.g. for a GW-BSE calculation in a polarizable
environment

.. code-block:: text

   Please cite the methods used in this run:
     [1] Baumeier, Wehner, Renaud, ..., VOTCA: multiscale frameworks ...
     [2] Ruehle et al., Microscopic simulations of charge transport ...
     [3] Thole, Chem. Phys. 59, 341 (1981)
     ...

``--cite <name>`` prints the references of a calculator or tool without
running it, including those that depend on its options
(e.g. ``xtp_tools --cite dftgwbse``). The references are

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - References
   * - VOTCA (every run)
     - [Baumeier:2024]_
   * - Structure-based calculators, KMC
     - [Ruehle:2011]_; exciton KMC (``kmclifetime``) [Wehner:2020]_
   * - Graph random walk (``grw``, ``latticegen``, KMC solvers ``exact``/``grw``)
     - [Chen:2025a]_, [Stenzel:2014]_
   * - Stochastic network model, traps
     - [Baumeier:2012]_, [Kordt:2014]_; traps [Chen:2025b]_
   * - DFT couplings, projection on dimer orbitals
     - [Baumeier:2010]_; POD [Ghan:2020]_
   * - GW-BSE
     - [Wehner:2018]_, [Tirimbo:2020a]_; QSGW [vanSchilfgaarde:2006]_,
       SRG regularization [Marie:2023]_
   * - BSE couplings, diabatization
     - [Wehner:2017]_; [Tirimbo:2024b]_
   * - DFT options
     - constrained DFT [Wu:2005]_; second-order SCF [Helmich-Paris:2021]_;
       projection-based embedding [Manby:2012]_, [Sundaram:2024]_
   * - Embedded GW-BSE (``qmmm``)
     - [Baumeier:2014]_, [Tirimbo:2020b]_, [Tirimbo:2024a]_; polar region
       [Thole:1981]_; Ewald region and ``ewaldbackground``
       [Poelking:2016]_, [Poelking:2015]_
   * - CHELPG partial charges
     - [Breneman:1990]_


Development
-----------

The core development is currently taking place at the Eindhoven University of Technology, Netherlands, and Los Alamos National Laboratory.


Copyright
---------

VOTCA is free software. The entire package is available under the Apache
License. For details, check the LICENSE file in the source code. The
VOTCA source code is available on `GitHub <https://github.com/votca/votca>`__.



