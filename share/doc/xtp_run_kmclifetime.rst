kmclifetime
***********
Calculator to perform the Kinetic Monte Carlo simulations of singlets with decay. The following table contains the default input options for the calculator,

.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - numberofinsertions
     - 
     - 
     - 
   * - seed
     - Integer to initialise the random number generator
     - 23
     - int+
   * - numberofcarriers
     - Number of electrons/holes in the simulation box
     - 1
     - int+
   * - injectionpattern
     - Name pattern that specifies on which sites injection is possible. Use the wildcard '*' to inject on any site.
     - *
     - 
   * - lifetimefile
     - File from which lifetimes are read in.
     - lifetimes.xml
     - 
   * - temperature
     - Temperature in Kelvin.
     - 300
     - int+
   * - trajectoryfile
     - Name of the trajectory file
     - 
     - 
   * - carrierenergy.run
     - Switch on/off
     - false
     - bool
   * - carrierenergy.outputfile
     - File to write energies to
     - energies.csv
     - 
   * - carrierenergy.alpha
     - Smoothing energy profile
     - 0.05
     - float+
   * - carrierenergy.outputsteps
     - Write to file every x steps
     - 10
     - int+


.. Note::
  An *xml* file containing the defaults for the `kmclifetime` calculator can be found at `${VOTCASHARE}/xtp/xml/kmclifetime.xml`.
