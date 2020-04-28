excitoncoupling
***************  
Exciton couplings from serialized orbital files.

The following table contains the input options for the calculator,

.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - classical
     - Exciton couplings from classical transition charges
     - 0
     - int+
   * - mpsA
     - classical transition charges for segment A
     - 
     - 
   * - mpsB
     - classical transition charges for segment B
     - 
     - 
   * - output
     - Output file
     - votca_excitoncoupling.xml
     - 
   * - bsecoupling_options.spin
     - Spin type for couplings
     - singlet
     - singlet,triplet,all
   * - bsecoupling_options.degeneracy
     - Criterium for the degeneracy of two levels
     - 0
     - 
   * - moleculeA.states
     - Number of excitons considered
     - 5
     - int+
   * - moleculeA.occLevels
     - occupied levels for CTstates
     - 5
     - int+
   * - moleculeA.unoccLevels
     - unoccupied levels for CTstates
     - 5
     - int+
   * - moleculeB.states
     - Number of excitons considered
     - 5
     - int+
   * - moleculeB.occLevels
     - occupied levels for CTstates
     - 5
     - int+
   * - moleculeB.unoccLevels
     - unoccupied levels for CTstates
     - 5
     - int+
   * - orbitalsA
     - Serialized orbitals file
     - A.orb
     - 
   * - orbitalsB
     - Serialized orbitals file
     - B.orb
     - 
   * - orbitalsAB
     - Serialized orbitals file
     - AB.orb
     -
     

.. note::
   An *xml* file containing the defaults for the `excitoncoupling` calculator can be found at `${VOTCASHARE}/xtp/xml/excitoncoupling.xml`.
