spectrum
********
Convolutes singlet spectrum with gaussian or lorentzian function.

The following table contains the input options for the calculator,

+----------------+------------------------------------+-------------------+--------------------+
| Property Name  |  Description                       | Default Value     |   Valid Input      |
+================+====================================+===================+====================+
|          fwhm  |          peak width in eV          |        0.2        |                    |
+----------------+------------------------------------+-------------------+--------------------+
|         lower  |   lower bound of spectrum in eV    |        0.0        |                    |
+----------------+------------------------------------+-------------------+--------------------+
|         upper  |   upper bound of spectrum in eV    |        3.5        |                    |
+----------------+------------------------------------+-------------------+--------------------+
|         points | | datapoints between upper and     |        100        |                    |
|                | | lower to calculate               |                   |                    |
+----------------+------------------------------------+-------------------+--------------------+
|          type  |print put energy/wavelength (eV/nm) |       enery       |                    |
+----------------+------------------------------------+-------------------+--------------------+
|         minexc | | lowest exciton to include in     |         0         |                    |
|                | | spectrum                         |                   |                    |
+----------------+------------------------------------+-------------------+--------------------+
|         maxexc | | highest exciton to include in    |       10000       |                    |
|                | | spectrum                         |                   |                    |
+----------------+------------------------------------+-------------------+--------------------+
|         shift  |   shift spectrum by amount of eV   |        0.0        |                    |
+----------------+------------------------------------+-------------------+--------------------+


.. note::
   An *xml* file containing the defaults for the `spectrum` calculator can be found at `${VOTCASHARE}/xtp/xml/spectrum.xml`.
