# About the Tutorials
This repository contains several tutorials for the **VOTCA-XTP** package
that show you how to perform the following calculations:

*  	A point energy calculation using DFT-GWBSE
*  	A geometry optimization using DFT-GWBSE
*  	A QMMM/DFT-GWBSE workflow using a GROMACS trajectory
*  	A QMMM/DFT-GWBSE workflow using a LAMMPS trajectory
*  	A Kinetic Monte-Carlo simulation for excited states using a LAMMPS trajectory
*  	A Kinetic Monte-Carlo simulation for excited states using a GROMACS trajectory

## Installation
*   In order to try then You will need to follow the **VOTCA** installation instructions
[here](https://github.com/votca/votca/blob/master/share/doc/INSTALL.md)
*   You will also need to install [Python and Jupyter](https://jupyter.readthedocs.io/en/latest/install.html)

Further information on VOTCA can be found at
http://www.votca.org

## Tutorials Location

### QMMM Tutorials
You can find the GROMACS based tutorial at:

``GROMACS/Methane/QMMM_GROMACS.ipynb``

and the one based on LAMMPS at:

``LAMMPS/Thiophene/QMMM_LAMMPS.ipynb``

### DFT-GWBSE TutorialS
The single point calculation and optmization tutorials are located at:

``
tools/dftgwbse_CH4/DFTGWBSE_ENERGY.ipynb
tools/dftgwbse_CO_geoopt/DFTGWBSE_OPTIMIZATION.ipynb
``

### KMC Tutorials
the KMC tutorials using GROMACS and LAMMPS trajectories, can be found at:

``
GROMACS/KMC_Methane/GROMACS_KMC.ipynb
LAMMPS/KMC_Thiophene/LAMMPS_KMC.ipynb
``

## Citation

The development of VOTCA is mainly funded by academic research grants. If you
use this package, please cite the VOTCA papers:

*   Excited-State Electronic Structure of Molecules Using Many-Body Green’s Functions: Quasiparticles and Electron-Hole Excitations with VOTCA-XTP,
    G. Tirimbo, V. Sundaram, O. Caylak, W. Scharpach, J. Sijen, C. Junghans, J. Brown, F. Zapata Ruiz, N. Renaud, J. Wehner, and B. Baumeier,
    ChemRxiv:11477895 (2020).

*   Electronic Excitations in Complex Molecular Environments: Many-Body Green’s
    Functions Theory in VOTCA-XTP Jens Wehner, Lothar Brombacher, Joshua Brown,
    Christoph Junghans, Onur Caylak, Yuriy Khalak, Pranav Madhikar, Gianluca
    Tirimbo, Björn Baumeier J. Chem. Theory Comput. 14, 6353 (2018).

*   Microscopic simulations of charge transport in disordered organic semiconductors
    V. Ruehle, A. Lukyanov, F. May, M. Schrader, T. Vehoff, J. Kirkpatrick, B. Baumeier and D. Andrienko
    J. Chem. Theo. Comp. 7, 3335-3345 (2011) 

*   Versatile Object-oriented Toolkit for Coarse-graining Applications,
    V.Ruehle, C. Junghans, A. Lukyanov, K. Kremer, D. Andrienko,
    J. Chem. Theo. Comp. 5 (12), 3211 (2009) 

In case of questions, please post them in the google discussion group
for votca at [here](https://groups.google.com/forum/#!forum/votca)

You can contact the VOTCA Development Team at devs@votca.org.



