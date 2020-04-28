coupling
********
Compute the electronic couplings from log and orbital files.
The following table contains the input options for the calculator,

.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - package.name
     - Name of the DFT package
     - xtp
     - xtp,orca
   * - dftcoupling_options.degeneracy
     - Criterium for the degeneracy of two levels
     - 0
     - float+
   * - dftcoupling_options.levA
     - Output HOMO, ..., HOMO-levels; LUMO, ..., LUMO+levels, molecule A
     - 1
     - int+
   * - dftcoupling_options.levB
     - Output HOMO, ..., HOMO-levels; LUMO, ..., LUMO+levels, molecule B
     - 1
     - int+
   * - moleculeA.log
     - Log file of molecule A
     - A.log
     - 
   * - moleculeA.orbitals
     - Orbitals file
     - A.orb
     - 
   * - moleculeB.log
     - Log file of molecule B
     - B.log
     - 
   * - moleculeB.orbitals
     - Orbitals file
     - B.orb
     - 
   * - dimerAB.log
     - Log file of dimer AB
     - AB.log
     - 
   * - dimerAB.orbitals
     - Orbitals file
     - A.orb
     - 

.. note::
   An *xml* file containing the defaults for the `coupling` calculator can be found at `${VOTCASHARE}/xtp/xml/coupling.xml`.
