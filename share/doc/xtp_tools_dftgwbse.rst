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
|    optimize         | Perform a Geometry optimization    |  false            | bool               |
+---------------------+------------------------------------+-------------------+--------------------+
|    read_guess       | Read Wave function guess           | false             | bool               |
+---------------------+------------------------------------+-------------------+--------------------+
| write_charges       | Print computed charges             | false             | bool               |
+---------------------+------------------------------------+-------------------+--------------------+

**GWBSE** options:

+------------------------------+------------------------------------+---------------------+--------------------+
|  Property Name               |  Description                       | Default Value       |   Valid Input      |
+==============================+====================================+=====================+====================+ 
|      gwbse_engine.tasks      |    Task to perform                 |input,dft,parse,gwbse|                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|         gwbse.ranges         | | default: all levels in RPA,      |  default            | | default, factor, |
|                              | | 1:2*HOMO in QP and all in BSE;   |                     | | explicit, full   |
|                              | | other options: factor,explicit   |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|         gwbse.rpamax         | | only needed, if ranges is factor |                     |                    |
|                              | | or explicit, number of levels in |                     |                    |
|                              | | rpa                              |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|         gwbse.qpmin          | | only needed, if ranges is factor |                     |                    | 
|                              | | or explicit, lowest MO to be     |                     |                    |
|                              | | used in GW                       |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|         gwbse.qpmax          | | only needed, if ranges is factor |                     |                    |
|                              | | or explicit, highest MO to be    |                     |                    |
|                              | | used in GW                       |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|         gwbse.bsemin         | | only needed, if ranges is factor |                     |                    |
|                              | | or explicit, lowest MO to be     |                     |                    |
|                              | | used in BSE                      |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|         gwbse.bsemax         | | only needed, if ranges is factor |                     |                    |
|                              | | or explicit, highest MO to be    |                     |                    |
|                              | | used in BSE                      |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|        gwbse.vxc.grid        | | grid quality                     |      medium         | | xcoarse, coarse, |
|                              |                                    |                     | | medium,fine,     |
|                              |                                    |                     | | xfine            |
+------------------------------+------------------------------------+---------------------+--------------------+
|     gwbse.scissor_shift      | | preshift unoccupied MOs by a     |        0.0          |                    |
|                              | | constant for GW calculation      |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|          gwbse.mode          | | use single short (G0W0) or       |       evGW          |   G0W0, evGW       |
|                              | | self-consistent GW (evGW)        |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|         gwbse.tasks          | | tasks to do gw,singlets,triplets |    gw,singlets      | | gw,singlets,     |
|                              | | or all                           |                     | | triplets or all  |
+------------------------------+------------------------------------+---------------------+--------------------+
|    gwbse.sigma_integrator    | | self-energy correlation          |        ppm          |   ppm, exact       |
|                              | | integration method: ppm, exact   |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|          gwbse.eta           | | small parameter eta of the       |       1e-3          |                    |
|                              | | Green's function                 |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|       gwbse.qp_solver        | | QP equation solve method:        |       grid          |  fixedpoint, grid  |
|                              | | fixedpoint, grid                 |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|     gwbse.qp_grid_steps      |      number of QP grid points      |       1001          |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|    gwbse.qp_grid_spacing     |  spacing of QP grid points in Ha   |       0.001         |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|        gwbse.exctotal        |  | maximum number of BSE states to |        30           |                    |
|                              |  | calculate                       |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|         gwbse.useTDA         |  use TDA for BSE default `false`   |         0           |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|   gwbse.ignore_corelevels    | | exclude core MO level from       |        no           |                    |
|                              | | calculation on RPA,GW or BSE     |                     |                    |
|                              | | level                            |                     |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|  gwbse.gw_sc_max_iterations  |Maximum number of iterations in gvGW|        50       	  |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|    gwbse.gw_mixing_order     | | Mixing of QP energies in evGW -  |        20       	  |                    |
|                              | | 0: plain, 1: linear, >1:         |                 	  |                    |
|                              | | Anderson                         |                 	  |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|  gwbse.g_sc_max_iterations   |                                    |        100      	  |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
|  gwbse.use_Hqp_offdiagonals  | | Using symmetrized off-diagonal   |       false     	  |                    |
|                              | | elements of QP Hamiltonian in    |                 	  |                    |
|                              | | BSE                              |                 	  |                    |
+------------------------------+------------------------------------+---------------------+--------------------+
| gwbse_engine.redirect_logger |                                    |         0       	  |                    |
+------------------------------+------------------------------------+---------------------+--------------------+

**GWBSE Eigensolver** options:

+------------------------------+------------------------------------+-----------------+--------------------+
|  Property Name               |  Description                       | Default Value   |   Valid Input      |
+==============================+====================================+=================+====================+ 
|    dodavidson                |        use davidson solver         |         1       |      bool          |
+------------------------------+------------------------------------+-----------------+--------------------+
| davidson_correction          |                                    |        DPR      | DPR or OHLSEN      |
+------------------------------+------------------------------------+-----------------+--------------------+
|davidson_tolerance            |        loose,normal,strict         |      strict     | | loose,normal,    |
|                              |                                    |                 | | strict           |
+------------------------------+------------------------------------+-----------------+--------------------+
|  davidson_ortho              | orthogonalisation routine GS or QR |        GS       |   GS, QR           |
+------------------------------+------------------------------------+-----------------+--------------------+
| davidson_update              |  | how large the search space can  |       safe      |  min, max, safe    |
|                              |  | become min, safe, max           |                 |                    |
+------------------------------+------------------------------------+-----------------+--------------------+
| davidson_maxiter             |           max iterations           |        50       |                    |
+------------------------------+------------------------------------+-----------------+--------------------+
|   domatrixfree               | | solve without explicitly setting |         0       |                    |
|                              | | up BSE matrix, (slower but a lot |                 |  bool              |
|                              | | less memory required             |                 |                    |
+------------------------------+------------------------------------+-----------------+--------------------+


**Geometry Optimization** options:

+-------------------------+------------------------------------+-------------------+----------------------+
|  Property Name     	  |  Description              	       | Default Value     |   Valid Input        |
+=========================+====================================+===================+======================+
| state              	  | state to optimize                  | s1                |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
| statetracker       	  | property to track the state        |                   |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
| statetracker.filter	  |                                    |oscillatorstrength | | chargetransfer,    |
|                    	  |                                    |                   | | density,           |
|                    	  |                                    |                   | | localisation,      |
|                    	  |                                    |                   | | oscillatorstrength,|
|                    	  |                                    |                   | | overlap            |
+-------------------------+------------------------------------+-------------------+----------------------+
|   convergence.energy    |       default: 1.e-6 Hartree       |       1.e-6       |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
|  convergence.RMSForce   |    default: 3.e-5 Hartree/Bohr     |       3.e-5       |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
|  convergence.MaxForce   |    default: 1.e-4 Hartree/Bohr     |       1.e-4       |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
|  convergence.RMSStep    |        default: 6.e-4 Bohr         |       6.e-4       |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
|  convergence.MaxStep    |        default: 1.e-3 Bohr         |       1.e-3       |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
|    optimizer.method     |                                    |     BFGS-TRM      |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
|    optimizer.trust      |  initial trustregion in Angstrom   |       0.01        |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
|     forces.method       | | finite differences method,       |      central      |                      |
|                         | | central or forward               |                   |                      |
+-------------------------+------------------------------------+-------------------+----------------------+
|forces.CoMforce_removal  |   Remove total force on molecule   |         1         |                      | 
+-------------------------+------------------------------------+-------------------+----------------------+
|  forces.displacement    |      default: 0.001 Angstrom       |       0.001       |                      |
+-------------------------+------------------------------------+-------------------+----------------------+


.. Note::
   * *The `basisset`, `auxbasisset` and `functional` to run the DFT and GWBSE calculcations are taken from the *DFTGWBSE* section.
   * An *xml* file containing the defaults for the `dftgwbse` calculator can be found at `${VOTCASHARE}/xtp/xml/dftgwbse.xml`.
 
