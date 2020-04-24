partialcharges
**************
Tool to derive partial charges from QM results stores in serialized file.

The following table contains the input options for the calculator,

+------------------------------+------------------------------------+-------------------+--------------------+
|  Property Name               |  Description                       | Default Value     |   Valid Input      |
+==============================+====================================+===================+====================+
|     esp2multipole.state      | | ground-,excited or               |       n2S1        |   | groundstate    |
|                              | | transitionstate                  |                   |   | n2s1           |
|                              | | (groundstate,singlet,n2s1)       |                   |   | singlet        |
+------------------------------+------------------------------------+-------------------+--------------------+
|     esp2multipole.method     | | Method to use derive partial     |      CHELPG       |   | CHELPG         |
|                              | | charges, CHELPG and Mulliken     |                   |   | Mulliken       |
|                              | | implented                        |                   |                    |
+------------------------------+------------------------------------+-------------------+--------------------+
|    esp2multipole.gridsize    | | Grid accuracy for numerical      |       fine        |    | coarse        |
|                              | | integration within CHELPG and    |                   |    | medium        |
|                              | | GDMA coarse,medium,fine          |                   |    | fine          |
+------------------------------+------------------------------------+-------------------+--------------------+

.. note::
   An *xml* file containing the defaults for the `partialcharges` calculator can be found at `${VOTCASHARE}/xtp/xml/partialcharges.xml`.
