********
The architecture of VOTCA
********

VOTCA consists of three main libraries **tools**, **csg** and **xtp**. While the structure has partially historic reasons, in general the following rule of thumb holds. 

1. **tools** provides low level routines and in general things that are used in **csg** and **xtp**, e.g. the :code:`calculator` class, which is the abstract base class for all calculators used. The most used class in **csg** and **xtp** is the :code:`Property` class, which is responsible for reading and writing xml files and handles most options handling inside VOTCA. In general all things that are used in **csg** and **xtp** but are not domain related should end up in tools. Additionally **tools** pulls in the `Eigen <https://eigen.tuxfamily.org>`_ library which provides all Vector/Matrix/Linear algebra functionality.

2. **csg** contains everything that is used for coarse-graining atomistic potentials. This includes fileparsers and writers for many molecular dynamics fileformats, from *xyz* to `gromacs <http://www.gromacs.org/>`_ and `lammps <https://www.lammps.org/>`_ input formats. Next to the atom coordinates and velocities it also provides topology readers, which allow you to read in which atom belongs to which molecule and which interactions exist, etc... This information is then contained in the :code:`Topology` class, which is also the backbone of most algorithms. Next to C++ code, **csg** contains a variety of additional scripts in various scripting languages.

3. **xtp** contains functionality for exciton transport simulations on molecular dynamics structures (it uses the functionality provided by **csg** here) and nowadays a lot of functionality for electronic state calculation using DFT and GW-BSE and their coupling to classical environments, i.e. QM/MM. Similar to **csg**, it contains file readers for Quantum Chemistry packages (`ORCA <https://orcaforum.kofo.mpg.de/app.php/portal>`_ and its own DFT implementation). In contrast to **csg** it is OPENMP parallelized and also uses CUDA acceleration for some parts.

Besides the 3 main folders, we also have the **csg-tutorials** and **xtp-tutorials** folder, which contain the respective tutorials. 
