mapchecker
**********
Outputs .pdb files for segments, qmmolecules and classical segments to check the mapping.
The following table contains the default input options for the calculator,

.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - qm_states
     - qmstates states to check
     - n
     - 
   * - mp_states
     - multipole states to check
     -  e h
     - 
   * - map_file
     - xml file with segment definition
     - votca_map.xml
     - 
   * - md_pdbfile
     - PDB file with the MD calculation
     - md_segments.pdb
     - 
   * - qm_pdbfile
     - PDB file with the QM calculation
     - qm_segments.pdb
     - 
   * - mp_pdbfile
     - PDB file with the MP segments
     - mp_segments.pdb
     - 

.. Note::
  An *xml* file containing the defaults for the `mapchecker` calculator can be found at `${VOTCASHARE}/xtp/xml/mapchecker.xml`.
