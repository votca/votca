Installation
============

General (Source) Installation Instructions
------------------------------------------

To install the full package:

::

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'stable' or 'v1.4.1'
    git clone -b ${version} --recursive https://github.com/votca/votca.git
    cmake -Bbuilddir -DBUILD_CSGAPPS=ON -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} votca
    cmake --build builddir --parallel <number of cores>
    cmake --build builddir --target install


Dependency Installation
~~~~~~~~~~~~~~~~~~~~~~~

`Ubuntu <https://github.com/votca/buildenv/blob/master/ubuntu#L10-L13>`_

`Fedora <https://github.com/votca/buildenv/blob/master/fedora#L10-L15>`_

If you have an older version of a distribution cmake can automatically install the correct Gromacs and libint version by
adding ``-DBUILD_OWN_GROMACS_OWN_GROMACS=ON`` or ``-DBUILD_OWN_GROMACS_OWN_LIBINT=ON`` to your cmake command. 

On Ubuntu 20.04. you will not get the libecpint package. So in this case run `this <https://github.com/votca/buildenv/blob/55c88641046dbb6152cf7d9e536e17f29205479f/ubuntu#L30-L33>`_.


Resolving the 'not found' dependency errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming all the `dependencies <#dependency-installation>`__ have been
correctly installed, one or more might still appear 'not found' upon
configuring with ``cmake`` command (see above). In this case you will
need to find out the 'non-standard' location for each missed dependency
(most often a shared or dynamically loaded library, e.g.
``libgromacs.so.*``, ``libhdf5.so.*`` etc).

Error messages produced by Cmake usually provide instructive suggestions
for resolving dependency issues. In particular, an appropriate extra
``-D`` flag is necessary to specify the path to a missed package. You
will have to rerun the ``cmake`` command with the relevant flag(s)
added. For example, in the case of a locally installed version of
Gromacs:

::

    cmake -DBUILD_CSGAPPS=ON -DCMAKE_INSTALL_PREFIX=${prefix} -DGROMACS_INCLUDE_DIR=$HOME/gromacs/include -DGROMACS_LIBRARY=$HOME/gromacs/lib/libgromacs.so ..

Be careful to use exactly the option suggested in the error message! You
can also add ``-LH`` or ``-LAH`` options to the ``cmake`` command in
order to see the available options with brief explanations (note that
*changing some of the variables may result in more variables being
created*; run ``man cmake`` for more info).

*Only for Linux*: For each dependency package not found by Cmake
initially, it might be necessary to add the location of its ``lib``
directory to the environment variable ``LD_LIBRARY_PATH``, **before**
building and installing VOTCA, i.e. before running any ``make`` command.
For example:

::

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/gromacs/lib:$HOME/anaconda/lib

Note that ``LD_LIBRARY_PATH`` also needs to be set every time when
running an executable from the VOTCA installation afterwards (which can
be automated via user's login profile, e.g. in .bashrc). Alternatively,
CMake has options to *remember* where libraries came from at link time,
which can be enabled by setting ``CMAKE_INSTALL_RPATH_USE_LINK_PATH`` to
``ON``. VOTCA has enabled this option and a couple of other rpath
related tweaks when setting ``ENABLE_RPATH_INJECT`` to ``ON``.

Common CMake Flags
~~~~~~~~~~~~~~~~~~

-  ``BUILD_CSGAPPS`` - Build the extra csg applications repo (ON/OFF,
   Default OFF)
-  ``BUILD_XTP`` - Build the xtp repo (ON/OFF, Default OFF)
-  ``CMAKE_INSTALL_PREFIX`` - where to install the votca executables
   (Default is /usr/local/bin)
-  ``ENABLE_TESTING`` - compile tests (ON/OFF, Default OFF)

Other CMake Flags
~~~~~~~~~~~~~~~~~

-  ``BUILD_CSG_MANUAL`` - Build csg pdf manual
-  ``CMAKE_DISABLE_FIND_PACKAGE_<name>`` - Disable using an optional
   package called ``<name>`` (ON/OFF)
-  ``CMAKE_DISABLE_FIND_PACKAGE_HDF5`` - Disable using the optional
   package ``HDF5`` (ON/OFF, Default OFF; relevant only for the
   ``master`` branch)

Packages for various Linux Distributions
----------------------------------------

Fedora
~~~~~~

::

    dnf install votca-csg votca-xtp

CentOs
~~~~~~

::

    yum install epel-release
    yum update
    yum install votca-csg

OpenSuse
~~~~~~~~

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

    emerge votca-csg votca-xtp

Spack
~~~~~

`Spack <https://spack.io/>`__ is an package manager, which allows to
build VOTCA and all its dependencies:

::

    git clone clone https://github.com/spack/spack.git spack
    source spack/share/spack/setup-env.sh
    spack install votca-csg
    spack install votca-xtp

Stable version
^^^^^^^^^^^^^^

Spack can also install the latest stable version from git using:

::

    spack install votca-csg@stable

Development version
^^^^^^^^^^^^^^^^^^^

Spack can also install the latest development version from git using:

::

    spack install votca-csg@develop

Docker
~~~~~~

Votca is also available through docker and can be accessed and run with
the following docker commands:

::

    docker pull votca/votca
    docker run -it votca/votca /bin/bash

Release version
^^^^^^^^^^^^^^^

Docker can also install the latest released version, e.g.::

::

    docker run -it votca/votca:v1.6 /bin/bash

Stable version
^^^^^^^^^^^^^^

Docker can also install the latest stable version from git using:

::

    docker run -it votca/votca:stable /bin/bash

FreeBSD
~~~~~~~

On FreeBSD votca can be installed from a binary package (recommended):

::

    pkg install votka

or it can be rebuilt and installed from the sources (slow):

::

    cd /usr/ports/science/votca && make install clean

Linking Error: Undefined reference to
-------------------------------------

This error can happen due to a multitude of reasons. You may have
forgotten to add paths to the ``LD_LIBRARY_PATH`` or forgotten to import
an environment module. In these cases deleting the ``CMakeCache.txt`` in
your ``build`` folder and rerunning ``cmake`` will help. Unfortunately
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
