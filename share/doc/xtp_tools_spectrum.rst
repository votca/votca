spectrum
********
Convolutes singlet spectrum with gaussian or lorentzian function.

The following table contains the input options for the calculator,

.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - fwhm
     - peak width in eV
     - 0.2
     - float+
   * - lower
     - lower bound of spectrum in eV
     - 0.0
     - float+
   * - upper
     - upper bound of spectrum in eV
     - 3.5
     - float+
   * - points
     - datapoints between upper and lower to calculate
     - 100
     - int+
   * - type
     - print put energy/wavelength (eV/nm)
     - enery
     - energy,wavelength
   * - minexc
     - lowest exciton to include in spectrum
     - 0
     - int+
   * - maxexc
     - highest exciton to include in spectrum
     - 10000
     - int+
   * - shift
     - shift spectrum by amount of eV
     - 0.0
     - float+


.. note::
   An *xml* file containing the defaults for the `spectrum` calculator can be found at `${VOTCASHARE}/xtp/xml/spectrum.xml`.
