Workflow and Calculator Structure
#################################

A typical workflow for a morphology from molecular dynamics (see the
QMMM_GROMACS tutorial) is

#. ``xtp_map -t topol.tpr -c conf.gro -s system.xml -f state.hdf5``: map the
   MD configuration to molecules (segments) with the mapping file
   ``system.xml``, which also refers to the QM geometries and the
   multipole/polarizable (MPS) files of each segment and state;
   ``xtp_run -e mapchecker`` writes the mapped structures as .pdb files for
   inspection.

   Other MD formats work the same way, only the input files differ. For
   LAMMPS, give the data file as topology and a dump file as trajectory, e.g.
   ``xtp_map -t system.data -c traj.dump -s system.xml -f state.hdf5 -i 99``
   for frame 99. The LAMMPS reader derives names that the mapping file has to
   use: atom types are the element (recognized from the mass) followed by
   the LAMMPS type number (``C63``, ``S23``), and molecules are named by
   their composition in alphabetical order (``C4H4S1`` for thiophene; a
   suffix ``-0``, ``-1``, ... distinguishes different molecules with the
   same composition). ``xtp_map -m`` writes a skeleton mapping file with
   these names for any input format.

#. ``xtp_run -e neighborlist``: pairs of neighbouring segments.
#. ``xtp_run -e einternal``: reads internal (gas-phase) site energies and
   reorganization energies (e.g. ``-c energies_file=system.xml``) into the
   state file.
#. Site energies in the environment, e.g. ``xtp_parallel -e qmmm`` with
   classical regions only, followed by ``xtp_run -e eanalyze``.
#. Electronic couplings: ``xtp_parallel -e eqm`` (single-molecule DFT/GW-BSE),
   then ``xtp_parallel -e iqm`` (couplings of all pairs), or
   ``xtp_parallel -e iexcitoncl`` for classical exciton couplings;
   ``xtp_run -e ianalyze`` summarizes them.
#. Transport: ``xtp_run -e kmcmultiple`` (bulk mobility), ``kmctof`` and
   ``grw`` (time of flight), ``kmclifetime`` (exciton diffusion), see
   :doc:`transport_simulations`.

The ``xtp_parallel`` calculators work in three steps: ``-j write`` writes a
job file with one job per segment or pair, ``-j run`` runs the jobs (``-t``
jobs at the same time, each with ``-x`` OpenMP threads; several processes can
work on the same job file, each taking ``-q`` jobs at a time), and ``-j read``
reads the results back into the state file.

.. toctree::
   xtp_map_overview
   xtp_parallel_overview
   xtp_run_overview
   xtp_tools_overview
