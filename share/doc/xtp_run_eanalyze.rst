eanalyze
********
Calculator to compute the histogram and correlation function of site energies and pair energy differences. The following table contains the default input options for the calculator,

.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - resolution_sites
     - Bin size for site energy histogram
     - 0.01
     - float+
   * - resolution_pairs
     - Bin size for pair energy histogram
     - 0.01
     - float+
   * - resolution_spatial
     - Bin size for site energy correlation
     - 0.01
     - float+
   * - states
     - states to analyze
     - e h s t
     - 

.. Note::
  An *xml* file containing the defaults for the `eanalyze` calculator can be found at `${VOTCASHARE}/xtp/xml/eanalyze.xml`.
