|Codacy Badge| |CI| |Docker| |DOI|

This is VOTCA, which provides the two following subpackages:

-  VOTCA-CSG, a library which provides tools to develop coarse-grained
   potentials from atomistic simulation data
-  VOTCA-XTP, a library providing the DFT+GW-BSE method to calculate 
   electronically excited states in single molecules, in molecular materials 
   with quantum-classical embedding, as well as electron, hole, and exciton 
   dynamics in atomistic MD-trajectories.

Basic installation (see performance advice below)
#################################################
::

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'stable' or 'v2021.2'
    git clone -b ${version} https://github.com/votca/votca.git
    cmake -B builddir -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} votca
    cmake --build builddir --parallel <number of cores>
    cmake --build builddir --target install

Using this code via docker
##########################
::

    docker pull votca/votca
    docker run -it votca/votca /bin/bash



Performance advice
##################
VOTCA-XTP relies in the Eigen library for vector-matrix operations, and a lot of 
performance can be gained by enabling vectorization and/or use of Intel's ``MKL`` 
as backend, which is automatically detected by ``CMake``. Below are some recommendations
for different architectures:

Intel Processors
****************
``g++``, ``clang``, and ``ipcx`` from the Intel OneAPI basekit give similar performance 
when used with the MKL. No special flags have to be supplied to ``CMake``.

If ``g++`` or ``clang`` are used without MKL, add ``-DCMAKE_CXX_FLAGS=-march=native`` to enable 
vectorization in ``Eigen``.

As a rough estimate, runtimes with vectorization and  ``gcc/clang`` are 30% shorter than without
vectorization. Use of ``MKL`` reduces them by another 50%. 

AMD Processors
**************
We recommend using ``g++`` or ``clang`` rather than an Intel compiler on AMD. Enable 
vectorization in ``Eigen`` with ``-DCMAKE_CXX_FLAGS=-march=native``. 

If you have Intel's MKL installed, and it is found by ``CMake`` make sure you 
use ``-DCMAKE_CXX_FLAGS=-march=core-avx2`` (provided the CPU supports ``AVX2``).

We observed about 50% shorter runtimes with vectorization enabled. Using the MKL 
can be slightly faster, but the gains (if there are any) are at best a few %.

CUDA support
************
If your system has a ``NVIDIA`` GPU, enable offloading of matrix operations 
by ``-DUSE_CUDA=ON``. 



More detailed information:
##########################

1. `Installation <share/doc/INSTALL.rst>`__
2. `Further Information <http://www.votca.org>`__
3. `Developers Guide <share/doc/DEVELOPERS_GUIDE.rst>`__
4. `VOTCA\_LANGUAGE\_GUIDE <share/doc/VOTCA_LANGUAGE_GUIDE.rst>`__
5. `Code of Conduct <share/doc/CODE_OF_CONDUCT.rst>`__

You find the VOTCA-CTP repository
`here <https://gitlab.mpcdf.mpg.de/votca/votca>`__

.. |Codacy Badge| image:: https://app.codacy.com/project/badge/Grade/b5567bfcf2c8411a8057c47fa7126781
   :target: https://www.codacy.com/gh/votca/votca?utm_source=github.com&utm_medium=referral&utm_content=votca/votca&utm_campaign=Badge_Grade
.. |CI| image:: https://github.com/votca/votca/workflows/CI/badge.svg?branch=master
   :target: https://github.com/votca/votca/actions?query=workflow%3ACI+branch%3Amaster
.. |Docker| image:: https://github.com/votca/votca/workflows/Docker/badge.svg?branch=master
   :target: https://github.com/votca/votca/actions?query=workflow%3ADocker+branch%3Amaster
.. |DOI| image:: https://zenodo.org/badge/75022030.svg
   :target: https://zenodo.org/badge/latestdoi/75022030
