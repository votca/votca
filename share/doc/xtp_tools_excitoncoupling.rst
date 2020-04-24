excitoncoupling
***************  
Exciton couplings from serialized orbital files.

The following table contains the input options for the calculator,

+------------------------------+------------------------------------+--------------------------+--------------------+
|  Property Name               |  Description                       | Default Value    	       |   Valid Input      |
+==============================+====================================+==========================+====================+
|          classical           | | Exciton couplings from classical |         0        	       |                    |
|                              | | transition charges               |                  	       |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|             mpsA             | | classical transition charges for |                  	       |                    |
|                              | | segment A                        |                  	       |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|             mpsB             | | classical transition charges for |                  	       |                    |
|                              | | segment B                        |                  	       |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|            output            |            Output file             |votca_excitoncoupling.xml |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|   bsecoupling_options.spin   | | Spin type for couplings,         |      singlet             |                    |
|                              | | singlet,triplet,all              |                          |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|bsecoupling_options.degeneracy| | Criterium for the degeneracy of  |         0                |                    |
|                              | | two levels                       |                          |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|       moleculeA.states       |   Number of excitons considered    |         5                |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|     moleculeA.occLevels      |    occupied levels for CTstates    |         5                |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|    moleculeA.unoccLevels     |   unoccupied levels for CTstates   |         5                |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|       moleculeB.states       |   Number of excitons considered    |         5                |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|     moleculeB.occLevels      |    occupied levels for CTstates    |         5                |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|    moleculeB.unoccLevels     |   unoccupied levels for CTstates   |         5                |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|          orbitalsA           |      Serialized orbitals file      |       A.orb              |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|          orbitalsB           |      Serialized orbitals file      |       B.orb              |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+
|          orbitalsAB          |      Serialized orbitals file      |      AB.orb              |                    |
+------------------------------+------------------------------------+--------------------------+--------------------+


.. note::
   An *xml* file containing the defaults for the `excitoncoupling` calculator can be found at `${VOTCASHARE}/xtp/xml/excitoncoupling.xml`.
