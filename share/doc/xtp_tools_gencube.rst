gencube
*******
Calculator to generate a **cube** file using the computed molecular orbitals.

The following table contains the input options for the calculator,

+---------------------+------------------------------------+-------------------+--------------------+
|  Property Name      |  Description                       | Default Value     |   Valid Input      |
+=====================+====================================+===================+====================+
|       padding       | How far the grid should start from |        6.5        |   postive real     |
|                     | the molecule                       |                   |                    |
+---------------------+------------------------------------+-------------------+--------------------+
|       xsteps        |     Gridpoints in x-direction      |        25         |   natural number   |
+---------------------+------------------------------------+-------------------+--------------------+
|       ysteps        |     Gridpoints in y-direction      |        25         |   natural number   |
+---------------------+------------------------------------+-------------------+--------------------+
|       zsteps        |     Gridpoints in z-direction      |        25         |   natural number   |
+---------------------+------------------------------------+-------------------+--------------------+
|        state        |  State to generate cube file for   |         N         |                    |
+---------------------+------------------------------------+-------------------+--------------------+
|       diff2gs       | | For excited states output        |                   |                    |
|                     | | difference to groundstate        |       false       |   bool             |
+---------------------+------------------------------------+-------------------+--------------------+
|        mode         | | new: generate new cube file,     |                   |                    |
|                     | | substract: substract to cube     |                   |                    |
|		      | | files specified below            |        new        |   new, substract   |
+---------------------+------------------------------------+-------------------+--------------------+


.. note::
   An *xml* file containing the defaults for the `gencube` calculator can be found at `${VOTCASHARE}/xtp/xml/gencube.xml`.
