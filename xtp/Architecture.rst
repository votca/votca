*****
Architecture of votca-xtp
*****

**xtp** has 4 major subfolders.

1. **include** contains nearly all the header files, these are also installed and can be included from other libraries
2. **scripts** contains a small number of python scripts, which do small operations on xml and hdf5 files. 
3. **share** contains information about `basissets <https://en.wikipedia.org/wiki/Basis_set_(chemistry)>`_ and `pseudopotentials <https://en.wikipedia.org/wiki/Pseudopotential>`_ which can be then used in calculations. They are saved in xml file format. Additionally **share** also contains xml files with the descriptions of all options for all calculators in **xtp** (If you do not know, what a calculator is, do the tutorial first).
4. **src**: **src/libxtp** contains all the source files for the header files in **include** and a bit more. **src/tests** contains all the unittests for **xtp** ( if you want to run them add :code:`-DENABLE_TESTING=ON` to your cmake config command go to the *builddir* and run :code:`ctest`. Lastly **src/tools** contains the executables, which you should have seen in the tutorial. 

For more information look at the developer guide on `the votca webpage <www.votca.org>`_ and maybe watch this `video <https://youtu.be/QPnbB0dQ7Pk>`_.
In the following I want to try to explain some major information flows/designs in **xtp**.

######
How electronic structure calculations are run


##### 
How options are parsed in xtp


#####
How statefiles are parsed in


#### 
How the cuda support works

