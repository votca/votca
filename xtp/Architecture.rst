*************************
Architecture of votca-xtp
*************************

**xtp** has 4 major subfolders.

1. **include** contains nearly all the header files, these are also installed and can be included from other libraries
2. **scripts** contains a small number of python scripts, which do small operations on xml and hdf5 files. 
3. **share** contains information about `basissets <https://en.wikipedia.org/wiki/Basis_set_(chemistry)>`_ and `pseudopotentials <https://en.wikipedia.org/wiki/Pseudopotential>`_ which can be then used in calculations. They are saved in xml file format. Additionally **share** also contains xml files with the descriptions of all options for all calculators in **xtp** (If you do not know, what a calculator is, do the tutorial first).
4. **src**: **src/libxtp** contains all the source files for the header files in **include** and a bit more. **src/tests** contains all the unittests for **xtp** (if you want to run them add :code:`-DENABLE_TESTING=ON` to your cmake config command go to the *builddir* and run :code:`ctest`. Lastly **src/tools** contains the executables, which you should have seen in the tutorial. 

For more information look at the developer guide on `the votca webpage <www.votca.org>`_ and maybe watch this `video <https://youtu.be/QPnbB0dQ7Pk>`_.
In the following I want to try to explain some major information flows/designs in **xtp**.

How electronic structure calculations are run
#############################################

So DFT calculations need as an input a basisset, which tells you which kind of basisfunctions belong to a certain element and we need a list of atoms, their element and coordinates. **xtp** comes with a number of basissets in the `share/xtp/basis_sets <https://github.com/votca/votca/tree/master/xtp/share/xtp/basis_sets>`_ directory. So the :code:`basisset` class loads the basisset, loading a list of atoms is normally done via the :code:`qmmolecule` class, which can read xyz files. Otherwise you can use a `QMMapper  <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/segmentmapper.h>`_ to convert a molecular dynamics molecule into a :code:`qmmolecule`. Combining a :code:`basisset` with a :code:`qmmolecule` gives you an `AOBasis <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/aobasis.h>`_, which contains the information about the shells/basisfunctions for each atom. For most functionality you also need an auxiliary basisset for the RI approximation, you find more information about that `here <https://research.tue.nl/en/publications/investigation-of-exciton-properties-in-organic-materials-via-many>`_. Using the :code:`AOBasis` you can create the `single electron matrices <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/aomatrix.h>`_ and `electron-electron integrals <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/ERIs.h>`_, which are mostly produced via the `libint <https://github.com/evaleev/libint>`_ and `libecpint <https://github.com/robashaw/libecpint>`_ libraries (the latter is for pseudopotentials). For the exchange correlation functionals we use the `libxc library <https://www.tddft.org/programs/libxc/>`_. To evaluate this you have to perform a numerical integration on a `grid <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/vxc_grid.h>`_ and `evaluate the functional <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/vxc_potential.h>`_ . For DFT you have to do this in a self-consisten way which is done in the `dftengine <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/dftengine.h>`_. 

If you instead do not want to use the internal DFT implementation but use ORCA, this is all wrapped in an abstract class `qmpackage <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/qmpackage.h>`_ , where you can choose ORCA or xtpdft, see for example `dftgwbse <https://github.com/votca/votca/blob/master/xtp/src/libxtp/tools/dftgwbse.cc>`_. The results of DFT calculations are stored in the `orbitals class <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/orbitals.h>`_

The GW-BSE calculations are wrapped in `gwbse <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/gwbse.h>`_. They typically consist of a GW and a BSE step.  The GW step is in the `GW class <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/gw.h>`_ which can choose different ways to evaluate the self-energy, Sigma, which is implemented via a `factory <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/sigmafactory.h>`_. The BSE is found in the `bse class <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/bse.h>`_. Technically the BSE is mostly an eigenvalue decomposition of a very large matrix, which we never explicitly build, but only calculate `its effects <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/bse_operator.h>`_ on a testvector in the `davidson solver <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/davidsonsolver.h>`_. The results of the GW and BSE calculations are stored in the `orbitals class <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/orbitals.h>`_

How options are parsed in xtp
#############################

So all calculators can be supplied options in two ways via the :code:`-o` flag and then giving an xml file with the options or via the :code:`-c` flag and then adding a string like :code:`TAG=VALUE` or a combination of both. In `xtp application <https://github.com/votca/votca/blob/master/xtp/src/libxtp/xtpapplication.cc>`_ first it is checkd if options are given via the  :code:`-o` flag and otherwise an empty :code:`Property` object is create, to which then the options from the :code:`-c` are added. Then a `Optionshandler <https://github.com/votca/votca/blob/master/tools/include/votca/tools/optionshandler.h>`_ is created with the location of the default option files. 
First the optionhandler loads the default option file for that calculator, e.g. `dftgwbse <https://github.com/votca/votca/blob/master/xtp/share/xtp/xml/dftgwbse.xml>`_ then it replaces al links in that file with the subpackages e.g. `dftpackage <https://github.com/votca/votca/blob/master/xtp/share/xtp/xml/subpackages/dftpackage.xml>`_ and then adds the useroptions. It checks if the useroptions exist  and if an errors are in there. Afterwards the options are forwarded to the specific calculator. 

How statefiles are parsed in
############################

Statefiles are only used in the **xtp_parallel** and **xrp_run** executables, which derive from `stateapplication <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/stateapplication.h>`_, which again derives from the above mentioned `xtp application <https://github.com/votca/votca/blob/master/xtp/src/libxtp/xtpapplication.cc>`. Their parsing is implemented in `statesaver <https://github.com/votca/votca/blob/master/xtp/src/libxtp/statesaver.cc>`_. Underlying all statefiles are the hdf5interfaces that are also used for the orbitals file and are implemented in `checkpointwriter <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/checkpointwriter.h>`_ and `checkopointreader <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/checkpointreader.h>`_.

How the cuda support works
##########################

CUDA support is optional and is enabled via :code:`-DUSE_CUDA=ON` in cmake and then enabled via `votca_xtp_config <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/votca_xtp_config.h.in>`_. CUDA is basically only used for matrix multiplication. The matrix class which wraps the `CUBLAS <https://docs.nvidia.com/cuda/cublas/index.html>`_ class is `CudaMatrix <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/cudamatrix.h>`_, which has ownwership of gpumemory. Calling the CUDA code inside openmp is handled via `OpenMP_CUDA class <https://github.com/votca/votca/blob/master/xtp/include/votca/xtp/openmp_cuda.h>`_, which also implements the specific commands for CPU and GPU calls. 

How code is executed
####################

In general you either call **xtp_parallel** **xtp_run** or **xtp_tools** which then call **xtpapplication** or **stateapplication** which creates the specific calculator gives it some options and runs it. 