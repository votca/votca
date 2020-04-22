dftgwbse
********

the `dftgwbse` calculator optionally takes an *xml* file that have the following hierarchical
structure:

::

   dftgwbse
   ├── general dftgwbse options
   ├── DFTpackage
   │   └── general dftpackage options
   ├── gwbse_engine
   │   ├── task
   │   └── gwbse options
   ├── geometry optimization
   │   └── geometry optimization options


An *xml* file containing the defaults for the `dftgwbse` calculator can be found at `${VOTCASHARE}/xtp/xml/dftgwbse.xml`.


General **DFTGWBSE** options:

+---------------------+------------------------------------+-------------------+--------------------+
|  Property Name      |  Description              	   | Default Value     |   Valid Input      |
+=====================+====================================+===================+====================+
|    mode             | Calculator mode           	   | energy            | energy, optimize   |
+---------------------+------------------------------------+-------------------+--------------------+
|    bassiset         | Basis set for MOs         	   | def2-tzvp         |                    |
+---------------------+------------------------------------+-------------------+--------------------+
|    auxbasisset      | Auxiliary basis set for RI	   | aux-def2-tzvp     |                    |
+---------------------+------------------------------------+-------------------+--------------------+
|    functional       | Functional name according to LIBXC | XC_HYB_GGA_XC_PBEH|                    |
+---------------------+------------------------------------+-------------------+--------------------+

**DFTPackage** options:

+---------------------+------------------------------------+-------------------+--------------------+
|  Property Name      |  Description              	   | Default Value     |   Valid Input      |
+=====================+====================================+===================+====================+
|     name            | Engine use to run DFT              | xtp               | xtp, orca          |
+---------------------+------------------------------------+-------------------+--------------------+
|     charge          | Molecular charge                   | 0                 | integer            |
+---------------------+------------------------------------+-------------------+--------------------+
|     spin            | Molecular multiplicity             | 1                 | natural number     |
+---------------------+------------------------------------+-------------------+--------------------+
|     read_guess      | Read Wave function guess           | false             | bool               |
+---------------------+------------------------------------+-------------------+--------------------+
| write_charges       | Print computed charges             | false             | bool               |
+---------------------+------------------------------------+-------------------+--------------------+

**GWBSE** options:

+---------------------+------------------------------------+-------------------+--------------------+
|  Property Name      |  Description              	   | Default Value     |   Valid Input      |
+=====================+====================================+===================+====================+
| ranges              |                                    |  default          |  default, factor,  |
|                     |                                    |                   |  explicit, full    |
+---------------------+------------------------------------+-------------------+--------------------+

**Geometry Optimization** options:

+---------------------+------------------------------------+-------------------+--------------------+
|  Property Name      |  Description              	   | Default Value     |   Valid Input      |
+=====================+====================================+===================+====================+
| state               | state to optimize                  | s1                |                    |
+---------------------+------------------------------------+-------------------+--------------------+
| statetracker        | property to track the state        |                   |                    |
+---------------------+------------------------------------+-------------------+--------------------+
| statetracker.filter |                                    |oscillatorstrength | chargetransfer,    |
|                     |                                    |                   | density,           |
|                     |                                    |                   | localisation,      |
|                     |                                    |                   | oscillatorstrength,|
|                     |                                    |                   | overlap            |
+---------------------+------------------------------------+-------------------+--------------------+


.. Note::
   The `basisset`, `auxbasisset` and `functional` to run the DFT and GWBSE calculcations are taken from the *DFTGWBSE* section.

