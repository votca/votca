partialcharges
**************
Tool to derive partial charges from QM results stores in serialized file.

The following table contains the input options for the calculator,

.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - esp2multipole.state
     - ground-,excited or transitionstate
     - n2S1
     - groundstate,singlet,n2s1
   * - esp2multipole.method
     - Method to use derive partial charges
     - CHELPG
     - CHELPG,Mulliken
   * - esp2multipole.gridsize
     - Grid accuracy for numerical integration within CHELPG and GDMA
     - fine
     - coarse,medium,fine


.. note::
   An *xml* file containing the defaults for the `partialcharges` calculator can be found at `${VOTCASHARE}/xtp/xml/partialcharges.xml`.
