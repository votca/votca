*****
Architecture of votca-xtp
*****

**xtp** has 4 major subfolders.

1. **include** contains nearly all the header files, these are also installed and can be included from other libraries
2. **scripts** contains a small number of python scripts, which do small operations on xml and hdf5 files. 
3. **share** contains information about `basissets <https://en.wikipedia.org/wiki/Basis_set_(chemistry)>`_ and `pseudopotentials <https://en.wikipedia.org/wiki/Pseudopotential>`_ which can be then used in calculations. They are saved in xml file format. Additionally **share** also contains xml files with the descriptions of all options for all calculators in **xtp** (If you do not know, what a calculator is, do the tutorial first).
4. **src**: **src/libxtp** contains all the source files for the header files in **include** and a bit more. **src/tests** contains all the unittests for **xtp** (if you want to run them add :code:`-DENABLE_TESTING=ON` to your cmake config command go to the *builddir* and run :code:`ctest`. Lastly **src/tools** contains the executables, which you should have seen in the tutorial. 

For more information look at the developer guide on `the votca webpage <www.votca.org>`_ and maybe watch this `video <https://youtu.be/QPnbB0dQ7Pk>`_.
In the following I want to try to explain some major information flows/designs in **xtp**.

######
How electronic structure calculations are run

So DFT calculations need as an input a basisset, which tells you which kind of basisfunctions belong to a certain element and we need a list of atoms, their element and coordinates. **xtp** comes with a number of basissets in the `share/xtp/basis_sets <https://github.com/votca/votca/tree/master/xtp/share/xtp/basis_sets>`_ direcotry. So the :code:`basisset`class loads the basisset, loading a list of atoms is normally done via the :code:`qmmolecule` class, which can read xyz files. Otherwise you can use a `QMMapper <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/segmentmapper.h>` to convert a molecular dynamics molecule into a :code:`qmmolecule`. Combining a :code:`basisset` with a :code:`qmmolecule` gives you an `AOBasis <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/aobasis.h>`, which contains the information about the shells/basisfunctions for each atom. For most functionality you also need an auxiliary basisset for the RI approximation, you find more information about that `here <https://research.tue.nl/en/publications/investigation-of-exciton-properties-in-organic-materials-via-many>`_. Using the :code:`AOBasis` you can create the `single electron matrices <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/aomatrix.h>`_ and `electron-electron integrals <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/ERIs.h>`_, which are mostly produced via the `libint <https://github.com/evaleev/libint>`_ and `libecpint <https://github.com/robashaw/libecpint>`_ libraries (the latter is for pseudopotentials). For the exchange correlation functionals we use the `libxc library <https://www.tddft.org/programs/libxc/>`_. To evaluate this you have to perform a numerical integration on a `grid <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/vxc_grid.h>`_ and `evaluate the functional <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/vxc_potential.h>`_ . For DFT you have to do this in a self-consisten way which is done in the `dftengine <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/dftengine.h>`_. 

If you instead do not want to use the internal DFT implementation but use ORCA, this is all wrapped in an abstract class `qmpackage <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/qmpackage.h>`_ , where you can choose ORCA or xtpdft, see for example `dftgwbse <https://github.com/votca/votca/blob/master/xtp/src/libxtp/tools/dftgwbse.cc>`_. The results of DFT calculations are stored in the `orbitals class <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/orbitals.h>`_

The GW-BSE calculations are wrapped in `gwbse https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/gwbse.h`_. They typically consist of a GW and a BSE step.  The GW step is in the `GW class <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/gw.h>`_ which can choose different ways to evaluate the self-energy, Sigma, which is implemented via a `factory <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/sigmafactory.h>`_. The BSE is found in the `bse class <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/bse.h>`_. Technically the BSE is mostly an eigenvalue decomposition of a very large matrix, which we never explicitly build, but only calculate `its effects <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/bse_operator.h>`_ on a testvector in the `davidson solver <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/davidsonsolver.h>`_. The results of the GW and BSE calculations are stored in the `orbitals class <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/orbitals.h>`_


##### 
How options are parsed in xtp


#####
How statefiles are parsed in


#### 
How the cuda support works

