ianalyze
********
Calculator to compute a squared logarithm histogram for the couplings. The following table contains the default input options for the calculator,

.. list-table:: Description
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - resolution_logJ2
     - Bin size of histogram log(J2)
     - 0.5
     - float+
   * - resolution_space
     - Bin size for r in log(J2(r))
     - 0.05
     - float+
   * - states
     - States for which to calculate the histogram
     - e h s t
     - 

.. Note::
  An *xml* file containing the defaults for the `ianalyze` calculator can be found at `${VOTCASHARE}/xtp/xml/ianalyze.xml`.
