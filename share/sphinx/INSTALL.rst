.. _install:

Installation
============

General (Source) Installation Instructions
------------------------------------------

To install the full package:

::

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'v2025'
    git clone -b ${version} https://github.com/votca/votca.git
    cmake -B builddir -S votca -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix}
    cmake --build builddir --parallel <number of cores>
    ctest --test-dir builddir
    cmake --build builddir --target install


Dependency Installation
~~~~~~~~~~~~~~~~~~~~~~~

`Ubuntu <https://github.com/votca/buildenv/blob/master/ubuntu#L10-L13>`_

`Fedora <https://github.com/votca/buildenv/blob/master/fedora#L10-L15>`_

If you have an older version of a distribution, cmake can be run with the ``-DBUILD_OWN_GROMACS=ON`` or ``-DBUILD_OWN_LIBINT=ON`` flags to automatically install the correct GROMACS and libint version.

On Ubuntu 20.10 and older you will not get version 1.0.5 of the libecpint package. So in this case run `this <https://github.com/votca/buildenv/blob/55c88641046dbb6152cf7d9e536e17f29205479f/ubuntu#L30-L33>`_ or
download the package from `here <http://ftp.us.debian.org/debian/pool/main/libe/libecpint/>`_.


Testing your build
~~~~~~~~~~~~~~~~~~

To test if your build actually works, you can run the VOTCA testsuite:

::

    ctest --test-dir <path-to-build>



Resolving the 'not found' dependency errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming all the `dependencies <#dependency-installation>`__ have been
correctly installed, one or more might still appear with the 'not found' status upon
configuring with the ``cmake`` command (see above). In this case, you will
need to find the 'non-standard' location for each missing dependency
(most often a shared or dynamically loaded library, e.g.
``libgromacs.so.*``, ``libhdf5.so.*`` etc).

Error messages produced by CMake usually provide instructive suggestions
for resolving dependency issues. In particular, an appropriate extra
``-D`` flag is necessary to specify the path to a missed package. You
will have to rerun the ``cmake`` command with the relevant flag(s)
added. For example, in the case of a locally installed version of
GROMACS:

::

    cmake -DCMAKE_INSTALL_PREFIX=${prefix} -DGROMACS_INCLUDE_DIR=$HOME/gromacs/include -DGROMACS_LIBRARY=$HOME/gromacs/lib/libgromacs.so -S ..

Please also ensure GROMACS was build with ``-DGMX_INSTALL_LEGACY_API=ON``.

Be careful to use exactly the option suggested in the error message! You
can also add ``-LH`` or ``-LAH`` options to the ``cmake`` command in
order to see the available options with brief explanations (note that
*changing some of the variables may result in more variables being
created*; run ``man cmake`` for more info).

*Only for Linux*: For each dependency package not found by CMake
initially, it might be necessary to add the location of its ``lib``
directory to the environment variable ``LD_LIBRARY_PATH``, **before**
building and installing VOTCA, i.e. before running any ``make`` command.
For example:

::

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/gromacs/lib:$HOME/anaconda/lib

Note that ``LD_LIBRARY_PATH`` also needs to be set every time when
running an executable from the VOTCA installation afterwards (which can
be automated via the user's login profile, e.g. in .bashrc). Alternatively,
CMake has options to *remember* where libraries came from at link time,
which can be enabled by setting ``CMAKE_INSTALL_RPATH_USE_LINK_PATH`` to
``ON``. VOTCA has enabled this option and a couple of other rpath
related tweaks when setting ``ENABLE_RPATH_INJECT`` to ``ON``.

Common CMake Flags
~~~~~~~~~~~~~~~~~~

-  ``INSTALL_CSGAPPS`` - Install the extra csg applications repo (ON/OFF,
   Default OFF)
-  ``BUILD_XTP`` - Build the xtp repo (ON/OFF, Default OFF)
-  ``CMAKE_INSTALL_PREFIX`` - where to install the votca executables
   (Default is /usr/local/bin)
-  ``BUILD_TESTING`` - compile tests (ON/OFF, Default ON)

Other CMake Flags
~~~~~~~~~~~~~~~~~

-  ``CMAKE_DISABLE_FIND_PACKAGE_<name>`` - Disable using an optional
   package called ``<name>`` (ON/OFF)
-  ``CMAKE_DISABLE_FIND_PACKAGE_HDF5`` - Disable using the optional
   package ``HDF5`` (ON/OFF, Default OFF; relevant only for the
   ``master`` branch)


Performance advice
~~~~~~~~~~~~~~~~~~
VOTCA-XTP relies on the Eigen library for vector-matrix operations, and a lot of 
performance can be gained by enabling vectorization and/or use of Intel's ``MKL`` 
as backend, which is automatically detected by ``CMake``. Below are some recommendations
for different architectures:

Intel Processors
^^^^^^^^^^^^^^^^
``g++``, ``clang``, and ``ipcx`` from the Intel OneAPI basekit give similar performance 
when used with the MKL. No special flags have to be supplied to ``CMake``.

If ``g++`` or ``clang`` are used, the compiler option ``-march=native`` is automatically injected
into the build. If you compile VOTCA on a heterogeneous cluster with different instruction sets,
this may cause the executables to not run. Override this by specifying 
``-DCMAKE_CXX_FLAGS=-mtune=native`` (at probably lower performance), or perform node-type 
specific builds.

As a rough estimate, runtimes with vectorization and  ``gcc/clang`` are 30% shorter than without
vectorization. Use of ``MKL`` reduces them by another 50%. 

AMD Processors
^^^^^^^^^^^^^^
We recommend using ``g++`` or ``clang`` rather than an Intel compiler on AMD. Vectorization 
in ``Eigen`` is automatically enabled by injection of ``-march=native``. See above comment
about heterogeneous envionments. 

If you have Intel's MKL installed, and it is found by ``CMake``, performance will be helped 
but by how much idepends on which architecture-specific implementation MKL picks at runtime. 
This can be affected by a `vendor lock-in <https://en.wikipedia.org/wiki/Math_Kernel_Library>`__, 
and work-arounds are documented `for Intel MKL on AMD Zen <https://danieldk.eu/Posts/2020-08-31-MKL-Zen.html>`__ 
and `in general <https://documentation.sigma2.no/jobs/mkl.html#using-mkl-efficiently>`__. 
We advise that you test this on your specific architecture.

CUDA support
^^^^^^^^^^^^
If your system has a ``NVIDIA`` GPU, enable offloading of matrix operations 
by ``-DUSE_CUDA=ON``. 

Packages for various Linux Distributions
----------------------------------------

Fedora
~~~~~~

::

    dnf install votca

Or in older versions of Fedora:

::

    dnf install votca-csg votca-xtp

CentOS
~~~~~~

::

    yum install epel-release
    yum update
    yum install votca

Or in older versions of CentOS:

::

    yum install votca-csg

openSUSE
~~~~~~~~

::

    zypper install votca

Or in older versions of openSUSE:

::

    zypper install votca-csg votca-xtp

SLES
~~~~

::

    SUSEConnect -p PackageHub/12.2/x86_64
    zypper install votca-csg

Debian / Ubuntu
~~~~~~~~~~~~~~~

::

    apt-get install votca-csg

Gentoo
~~~~~~

::

    emerge votca

Spack
~~~~~

`Spack <https://spack.io/>`__ is a package manager, which has
the capability of building VOTCA and all its dependencies:

::

    git clone https://github.com/spack/spack.git spack
    source spack/share/spack/setup-env.sh
    spack install votca

Development version
^^^^^^^^^^^^^^^^^^^

Spack can also install the latest development version from git using:

::

    spack install votca@master

Other build options
^^^^^^^^^^^^^^^^^^^

Spack has other options:

::

    spack info votca
    
One useful option is to build votca without xtp:

::

    spack install votca~xtp

Conda-forge
~~~~~~~~~~~

::

    conda install -c conda-forge votca

Docker
~~~~~~

VOTCA is also available through docker and can be accessed and run with
the following docker commands:

::

    docker pull votca/votca
    docker run -it votca/votca /bin/bash

Release version
^^^^^^^^^^^^^^^

Docker can also install the latest released version, e.g.:

::

    docker run -it votca/votca:v2021.2 /bin/bash

FreeBSD
~~~~~~~

On FreeBSD VOTCA can be installed from a binary package (recommended):

::

    pkg install votca

or it can be rebuilt and installed from the sources (slow):

::

    cd /usr/ports/science/votca && make install clean

Linking Error: Undefined reference to
-------------------------------------

This error can occur for a multitude of reasons. You may have
forgotten to add paths to the ``LD_LIBRARY_PATH`` or forgotten to import
an environment module. In these cases, deleting the ``CMakeCache.txt`` in
your ``build`` folder and rerunning ``cmake`` will help. Unfortunately,
another cause might be ABI incompability between libraries due to
different libraries being compiled with different compilers or compiler
versions. Click `here <https://github.com/ICRAR/shark/issues/1>`__ for
an example.

GCC and MKL : undefined symbol:
-------------------------------

This can happen with some GCC versions. Adding the
``-Wl,--no-as-needed`` to ``CMAKE_EXE_LINKER_FLAGS:STRING=`` to the
``CMakeCache.txt`` in the ``build`` directory can fix this. For more
information look
`here <https://software.intel.com/en-us/articles/symbol-lookup-error-when-linking-intel-mkl-with-gcc-on-ubuntu>`__
