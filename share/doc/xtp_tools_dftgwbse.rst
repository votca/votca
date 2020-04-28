dftgwbse
********
Calculator to compute lectronic Excitations using the **GW-BSE** approximation.

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


.. list-table:: General **DFTGWBSE** options
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input
   * - mode
     - Calculator mode
     - energy
     - energy,optimize
   * - basisset
     - Basis set for MOs
     - def2-tzvp
     - 
   * - auxbasisset
     - Auxiliary basis set for RI
     - aux-def2-tzvp
     - 
   * - functional
     - Functional name(s) according to LIBXC
     - XC_HYB_GGA_XC_PBEH
     - 


.. list-table:: General **DFTPackage** options
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - package.name
     - Name of the DFT package
     - xtp
     - xtp,orca
   * - charge
     - Molecular charge
     - 0
     - int+
   * - spin
     - Molecular multiplicity
     - 1
     - int+
   * - optimize
     - Perform a Geometry optimization
     - false
     - bool
   * - read_guess
     - Read Wave function guess
     - false
     - bool
   * - write_charges
     - Print computed charges
     - false
     - bool


.. list-table:: General **GWBSE** options
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - gwbse_engine.tasks
     - task to compute
     - input,dft,parse,gwbse
     - 
   * - gwbse.ranges
     - default: all levels in RPA, 1:2*HOMO in QP and all in BSE; other options: factor,explicit
     - 
     - default,factor,explicit
   * - gwbse.rpamax
     - only needed, if ranges is factor or explicit, number of levels in rpa
     - 
     - 
   * - gwbse.qpmin
     - only needed, if ranges is factor or explicit, lowest MO to be used in GW
     - 
     - 
   * - gwbse.qpmax
     - only needed, if ranges is factor or explicit, highest MO to be used in GW
     - 
     - 
   * - gwbse.bsemin
     - only needed, if ranges is factor or explicit, lowest MO to be used in BSE
     - 
     - 
   * - gwbse.bsemax
     - only needed, if ranges is factor or explicit, highest MO to be used in BSE
     - 
     - 
   * - vxc.grid
     - grid quality
     - medium
     - xcoarse,coarse,medium,fine,xfine
   * - gwbse.scissor_shift
     - preshift unoccupied MOs by a constant for GW calculation
     - 0.0
     - float
   * - gwbse.mode
     - use single short (G0W0) or self-consistent GW (evGW)
     - evGW
     - evGW,G0W0
   * - gwbse.tasks
     - tasks to do
     - gw,singlets
     - 
   * - gwbse.sigma_integrator
     - self-energy correlation integration method
     - ppm
     - ppm, exact
   * - gwbse.eta
     - small parameter eta of the Green's function
     - 1e-3
     - float+
   * - gwbse.qp_solver
     - QP equation solve method
     - grid
     - fixedpoint,grid
   * - gwbse.qp_grid_steps
     - number of QP grid points
     - 1001
     - int+
   * - gwbse.qp_grid_spacing
     - spacing of QP grid points in Ha
     - 0.001
     - float+
   * - gwbse.exctotal
     - maximum number of BSE states to calculate
     - 30
     - int+
   * - gwbse.useTDA
     - use TDA for BSE
     - false
     - bool
   * - gwbse.ignore_corelevels
     - exclude core MO level from calculation on RPA,GW or BSE level
     - no
     - yes,no
   * - gwbse.gw_sc_max_iterations
     - Maximum number of iterations in gvGW
     - 50
     - int+
   * - gwbse.gw_mixing_order
     - Mixing of QP energies in evGW - 0: plain, 1: linear, >1: Anderson
     - 20
     - int+
   * - gwbse.g_sc_max_iterations
     - What is this again?
     - 100
     - int+
   * - gwbse.use_Hqp_offdiagonals
     - Using symmetrized off-diagonal elements of QP Hamiltonian in BSE
     - false
     - bool
   * - gwbse_engine.redirect_logger
     - Redirect Logger
     - false
     - bool

.. list-table:: **GWBSE Eigensolver** options
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - eigensolver.dodavidson
     - use davidson solver
     - true
     - bool
   * - eigensolver.davidson_correction
     - Davidson correction method
     - DPR
     - DPR,OHLSEN
   * - eigensolver.davidson_tolerance
     - Numerical tolerance
     - strict
     - loose,normal,strict
   * - eigensolver.davidson_ortho
     - orthogonalisation routine: Gram–Schmidt or QR
     - GS
     - GS,QR
   * - eigensolver.davidson_update
     -  how large the search space
     - safe
     - min,safe,max
   * - eigensolver.davidson_maxiter
     - max iterations
     - 50
     - int+
   * - eigensolver.domatrixfree
     - solve without explicitly setting up BSE matrix, (slower but a lot less memory required
     - false
     - bool
 


.. list-table:: General **Geometry Optimization** options
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - geometry_optimization.state
     - initial state to optimize for
     - s1
     - 
   * - statetracker.oscillatorstrength
     - 
     - 0.0001
     - float+
   * - statetracker.filters
     - property to track
     - oscillatorstrength
     - chargetransfer,density,localisation,oscillatorstrength,overlap
   * - convergence.energy
     - default: 1.e-6 Hartree
     - 1.e-6
     - float+
   * - convergence.RMSForce
     - default: 3.e-5 Hartree/Bohr
     - 3.e-5
     - float+
   * - convergence.MaxForce
     - default: 1.e-4 Hartree/Bohr
     - 1.e-4
     - float+
   * - convergence.RMSStep
     - default: 6.e-4 Bohr
     - 6.e-4
     - float+
   * - convergence.MaxStep
     - default: 1.e-3 Bohr
     - 1.e-3
     - float+
   * - optimizer.method
     - 
     - BFGS-TRM
     - 
   * - optimizer.trust
     - initial trustregion in Angstrom
     - 0.01
     - float+
   * - forces.method
     - finite differences method, central or forward
     - central
     - 
   * - forces.CoMforce_removal
     - Remove total force on molecule
     - true
     - bool
   * - forces.displacement
     - default: 0.001 Angstrom
     - 0.001
     - float+	   


.. Note::
   * *The `basisset`, `auxbasisset` and `functional` to run the DFT and GWBSE calculcations are taken from the *DFTGWBSE* section.
   * An *xml* file containing the defaults for the `dftgwbse` calculator can be found at `${VOTCASHARE}/xtp/xml/dftgwbse.xml`.
 
