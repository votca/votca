XTP Tutorials
=============

The tutorials are Jupyter notebooks in ``xtp-tutorials`` (installed to
``share/votca/xtp-tutorials``):

-  ``DFTGWBSE_ENERGY``: DFT and GW-BSE for a single molecule with ``xtp_tools``.
-  ``single_point`` and ``geometry_optimization``: the same with the Python
   (ASE) interface; they are built only if the Python bindings are.
-  Three tutorials on the same morphology, an amorphous film of 1000
   thiophene molecules from LAMMPS (``thiophene``):

   -  ``MAPPING``: from the MD snapshot to the state file (mapping file,
      ``xtp_map``, ``mapchecker``, neighbour list).
   -  ``CHARGE_TRANSPORT``: reorganization energies, site energies in the
      periodic polarizable environment and couplings, then mobilities and
      times of flight.
   -  ``EXCITATIONS``: GW-BSE excitations of one molecule in vacuum, in a
      static and in a polarizable environment (iterative and embedded
      schemes).

-  ``TRANSPORT`` and ``STOCHASTIC_NETWORK``: the transport calculators on
   model systems (KMC, time of flight, graph random walks, stochastic
   network models).

.. toctree::
   DFTGWBSE_ENERGY
   geometry_optimization
   single_point
   MAPPING
   CHARGE_TRANSPORT
   EXCITATIONS
   TRANSPORT
   STOCHASTIC_NETWORK

   
