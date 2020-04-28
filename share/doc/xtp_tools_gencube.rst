gencube
*******
Calculator to generate a **cube** file using the computed molecular orbitals.

The following table contains the input options for the calculator,


.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - padding
     - How far the grid should start from the molecule
     - 6.5
     - float+
   * - xsteps
     - Gridpoints in x-direction
     - 25
     - int+
   * - ysteps
     - Gridpoints in y-direction
     - 25
     - int+
   * - zsteps
     - Gridpoints in z-direction
     - 25
     - int+
   * - state
     - State to generate cube file for
     - N
     - N
   * - diff2gs
     - For excited states output difference to groundstate
     - false
     - bool
   * - mode
     - new: generate new cube file, substract: substract to cube files specified below
     - new
     - new,substract


.. note::
   An *xml* file containing the defaults for the `gencube` calculator can be found at `${VOTCASHARE}/xtp/xml/gencube.xml`.
