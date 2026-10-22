.. _install:

Installation
============

General (Source) Installation Instructions
------------------------------------------

To install the full package:

::

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'v2026'
    git clone -b ${version} https://github.com/votca/votca.git
    cmake -B builddir -S votca -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix}
    cmake --build builddir --parallel <number of cores>
    ctest --test-dir builddir
    cmake --build builddir --target install


Dependency Installation
~~~~~~~~~~~~~~~~~~~~~~~

`Ubuntu <https://github.com/votca/buildenv/blob/master/ubuntu#L10-L13>`_

`Fedora <https://github.com/votca/buildenv/blob/master/fedora#L10-L15>`_

If you have an older version of a distribution, cmake can be run with the ``-DBUILD_OWN_LIBINT=ON`` flag to automatically install the correct libint version.

On Ubuntu 20.10 and older you will not get version 1.0.5 of the libecpint package. So in this case run `this <https://github.com/votca/buildenv/blob/55c88641046dbb6152cf7d9e536e17f29205479f/ubuntu#L30-L33>`_ or
download the package from `here <http://ftp.us.debian.org/debian/pool/main/libe/libecpint/>`_.

macOS (Homebrew)
^^^^^^^^^^^^^^^^

On macOS (Apple silicon, macOS 13.3 or newer) use Apple's ``clang`` from the
Xcode command line tools and the dependencies from `Homebrew <https://brew.sh>`__:

::

    xcode-select --install
    brew install cmake eigen boost libint libecpint libxc hdf5 fftw libomp gromacs open-babel pybind11
    cmake -B builddir -S votca -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} \
          -DCMAKE_PREFIX_PATH="$(brew --prefix libint)"
    cmake --build builddir --parallel $(sysctl -n hw.ncpu)

Notes:

-  Linear algebra uses Apple's Accelerate framework, which ``CMake`` finds
   by itself (``-- Found APPLE_ACCELERATE``). Do **not** install Homebrew's
   ``lapack`` or ``lapacke`` for VOTCA: they are not needed.
-  Apple's ``clang`` has no OpenMP runtime of its own. ``CMake`` takes
   Homebrew's ``libomp`` automatically (``-- Using Homebrew libomp: ...``),
   so ``OpenMP_ROOT`` does not need to be set. To use another OpenMP
   runtime, set ``OpenMP_ROOT`` (as a CMake or environment variable) to its
   prefix. Without any OpenMP runtime VOTCA builds, but runs single-threaded.
-  Python packages used by some tests and tools (``lxml h5py rdkit scipy
   pytest ase xmltodict``) are best installed in a virtual environment made
   with Homebrew's ``python3``.

See `Apple Silicon`_ below for how to check that Accelerate is used.


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
``g++``, ``clang``, and ``icpx`` from the Intel OneAPI basekit give similar performance 
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
about heterogeneous environments.

Intel's MKL also works on AMD processors and is found by ``CMake`` in the same way. Every DFT
and GW-BSE log names the backend; a build with MKL prints

::

    Using MKL overload for Eigen (BLAS and LAPACK)

The largest gain is in the symmetric eigensolvers and the multithreaded rank-k updates, which
Eigen's own implementation runs on a single core. On an 8-core AMD EPYC 7262 (Zen 2) with
oneMKL 2026.2, MKL's ``dsyevd`` diagonalized a 3000 x 3000 matrix in 1.1 s, against 25 s for
Eigen's own solver.

At runtime, MKL checks whether it runs on an Intel processor and otherwise takes slower code
paths. This check can be overridden by preloading a small library. Save as ``mkl_cpu_trick.c``

.. code-block:: c

    int mkl_serv_intel_cpu_true(void) { return 1; }

    /* MKL 2025.0 and later call this function instead */
    typedef int (*fakeintel_fptr)(void);
    fakeintel_fptr mkl_serv_get_cpu_true(void) { return &mkl_serv_intel_cpu_true; }

then compile it and preload it when running VOTCA:

::

    gcc -O2 -shared -fPIC -o libmkl_cpu_trick.so mkl_cpu_trick.c
    LD_PRELOAD=/path/to/libmkl_cpu_trick.so xtp_tools ...

On the system above, this made MKL's matrix products 18-48% faster and its eigensolvers
1.5 times faster. A naphthalene def2-TZVP evGW-BSE calculation took 5% less time
(DFT 11%, BSE 20%, GW 2%) with unchanged results. Older recipes that define only
``mkl_serv_intel_cpu_true``, or that set ``MKL_DEBUG_CPU_TYPE=5``, have no effect with current
MKL versions. The override relies on MKL internals that Intel does not document or support,
so it may stop working with a future release; check the MKL license terms before using it.
Background: `Intel MKL on AMD Zen <https://danieldk.eu/software/misc/intel-mkl-on-amd-zen>`__
and `Using MKL efficiently <https://documentation.sigma2.no/jobs/mkl.html#using-mkl-efficiently>`__.
We have only measured Zen 2, so we advise that you test this on your specific architecture.

AMD's own AOCL libraries are currently not supported as a backend by VOTCA. In the same
benchmark, they matched the patched MKL for matrix products but were 4-11 times slower in the
eigensolvers.

Apple Silicon
^^^^^^^^^^^^^
On macOS, Eigen's matrix products go to Apple's Accelerate framework, and the
symmetric eigensolvers and Cholesky factorizations call Accelerate's LAPACK
directly. Accelerate reaches the CPU's matrix unit, which Eigen's own kernels
(NEON, two doubles wide) cannot use. No flags are needed; ``-march=native``
(or ``-mcpu=native`` if the compiler rejects it) is injected as on Linux.

Every DFT and GW-BSE log names the backend. A correct build prints

::

    Using Apple Accelerate BLAS overload for Eigen, Accelerate LAPACK for eigensolvers and Cholesky

If the log says ``Using native Eigen implementation, no BLAS overload``,
Accelerate was not found, and large GW-BSE runs are several times slower:
delete the build directory and configure again.

Use all cores (the default, or ``-t`` equal to ``sysctl -n hw.ncpu``).
Limiting ``xtp`` to the performance cores did not make it faster on M1 and
M1 Pro machines: the efficiency cores still add to the threaded loops, and
Accelerate's speed does not depend on the thread count. Do not set
``VECLIB_MAXIMUM_THREADS``.

Apple GPUs cannot be used: they have no double-precision arithmetic.

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

MKL and OpenMP thread pinning
-----------------------------

With ``g++``, VOTCA links MKL's GNU threading layer (``MKL_THREADING=gnu_thread``), so only
GCC's OpenMP runtime ``libgomp`` is loaded. If MKL's Intel threading layer is used with ``g++``
instead (``-DMKL_THREADING=intel_thread``, or a build directory configured with an older VOTCA
version), Intel's OpenMP runtime ``libiomp5`` is loaded next to ``libgomp``. This works, but do
not set ``OMP_PROC_BIND`` or ``OMP_PLACES`` with such a build: ``libgomp`` then pins the main
thread to one core, and all of MKL's threads end up on that core (visible in ``htop`` as one
core at 100% and the rest idle). ``ldd $(which xtp_tools) | grep iomp5`` shows whether
``libiomp5`` is linked; configure a new build directory to switch.

GCC and MKL : undefined symbol:
-------------------------------

This can happen with some GCC versions. Adding the
``-Wl,--no-as-needed`` to ``CMAKE_EXE_LINKER_FLAGS:STRING=`` to the
``CMakeCache.txt`` in the ``build`` directory can fix this. For more
information look
`here <https://software.intel.com/en-us/articles/symbol-lookup-error-when-linking-intel-mkl-with-gcc-on-ubuntu>`__
