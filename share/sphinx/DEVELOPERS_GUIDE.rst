Developer and Contributor Guide
===============================

This page is designed to give new developers general guidelines for
implementing code consistent with the VOTCA and cpp style and standard.

-  `Reporting Bugs <#reporting-bugs>`__
-  `Making a release <#making-a-release>`__
-  `CPP Resources <#cpp-resources>`__
-  `CPP Coding Rules <#cpp-coding-rules>`__
-  `Testing <#testing>`__
-  `CPP Coding Style Guide <#cpp-coding-style-guide>`__
-  `CPP Comment Guide <#cpp-comment-guide>`__
-  `Failed Release Builds <#failed-release-builds>`__

Reporting Bugs
--------------

To report a bug, please create an issue on the appropriate GitHub repo.
Please be sure to provide as much information as possible such as:

-  The error messages
-  The operating system
-  What compiler was used
-  What dependencies were installed
-  The calculation that was being run

Issues can be directly created on the `GitHub repo <https://github.com/votca/votca/issues>`__.

Formatting code
---------------

VOTCA uses ``clang-format`` to format the code, code that is not
properly formatted is automatically rejected. The style file ``.clang-format`` is in the
root of the repository. CMake provides a ``format`` target (it currently covers only
``*.c`` and ``*.h`` files).
The easiest way to format your code is just a ``@votca-bot format`` comment in the PR, which then will automatically format your code.

Doxygen documentation
---------------------
A complete overview of all C++ classes and code can be found on https://doc.votca.org/.

VOTCA dev-tools
---------------

Running clang-format on every commit can be a drag, as can changing the
copyright in every header. Fortunately, you will find small scripts in the
`dev-tools repo <https://github.com/votca/dev-tools>`__, which can
automate this.

VOTCA Continuous Integration (GitHub Actions)
---------------------------------------------

Each pull request to master in the votca repository 
is built on a machine in the cloud using `GitHub actions <https://docs.github.com/en/actions>`__ 

VOTCA is built on various linux distributions, which are not all natively supported by GitHub actions. The builds therefore run inside
prebuilt `docker containers <https://www.docker.com/resources/what-container>`__ (``ghcr.io/votca/buildenv/<distribution>``), which contain all
dependencies of VOTCA (see :code:`buildenv` below).

The docker images can be found at `Docker Hub <https://hub.docker.com/u/votca>`__. The **votca/buildenv** containers are the basic containers, which contain all the dependencies VOTCA requires; VOTCA code itself is not included. The **votca/buildenv** can be found on `VOTCA's GitHub Container registry <https://github.com/orgs/votca/packages>`__. 
The actual containers used for running the test builds are built on top of the **votca/buildenv** containers, the resulting **votca/votca** container can be found on `Docker Hub <https://hub.docker.com/u/votca>`__ as well as `VOTCA's GitHub Container registry <https://github.com/orgs/votca/packages>`__.

More information can be found in the `GitHub workflow files <https://github.com/votca/votca/tree/master/.github/workflows>`__.

Making a Release
----------------

Similar to the VOTCA containers, releases are also handled by GitHub actions. :code:`votca/votca` has a :code:`release` workflow that can only be triggered manually.
To trigger it go `this GitHub Action <https://github.com/votca/votca/actions?query=workflow%3Arelease>`_. The release can only be made from the 
:code:`master` branch, but testing the creation of a release can be triggered on any branch. To make a release, trigger the action from the
:code:`master` branch, pick a new release tag in the :code:`release tag` box and type :code:`yesyesyes` into the deploy box. :code:`CHANGELOG.rst` must contain a section for the current development version (e.g. ``Version 2026-dev``) with at least one item; the release script renames it to the new version with the release date. A new release will trigger the creation of the release tag.

Release names
~~~~~~~~~~~~~

Some releases have names, so far we have:

-  1.1: SuperAnn - named after the spouse of a core developer
-  1.2: SuperDoris - named after the administrator at MPI-P (VOTCA's birthplace)
-  1.3: SuperUzma - named after the spouse of a core developer
-  1.4: SuperKurt - in occasion of Kurt Kremer's 60th birthday
-  1.5: SuperVictor - named after Victor Rühle, one of the original core developers
-  1.6: SuperPelagia - named after the spouse of a core developer
-  1.6.2: SuperGitta - in memory of the grandmother of a core developer


CPP Resources
-------------

A good starting point, is to take a look at the cpp standard. Though the
code has not always consistently followed the cpp standard we now make
an effort to really enforce it and follow best practices.

-  `Best
   Practices1 <https://github.com/cpp-best-practices/cppbestpractices>`__
-  `Best
   Practices2 <https://google.github.io/styleguide/cppguide.html>`__

CPP Coding Rules
----------------


Here are a few general rules that should be followed:

Files
~~~~~

-  Each class goes into a separate file.
-  Each filename should be the the name of the class it contains written in lowercase.

Includes
~~~~~~~~

-  When including a header file from within the same repo that you are
   working use the relative includes. This consists of using quotation
   marks i.e.

   #include "molecule.h"

-  When including from another module, for instance you are working
   in the csg module and want to include a file from the tools repo
   use the anglular brackets i.e.

   #include <votca/tools/molecule.h>

Header Files
~~~~~~~~~~~~

-  One class, one header.
-  When creating header guards use the template: VOTCA\_VOTCA-REPO-NAME\_CLASS-NAME\_H. Where
   "VOTCA-REPO-NAME" is replaced by whichever repo the header file is in, this could be
   tools, csg or xtp. The "CLASS-NAME" component should also be replaced, but by the name of the
   class described in the header file:

   #ifndef VOTCA\_VOTCA-REPO-NAME\_CLASS-NAME\_H #define
   VOTCA\_VOTCA-REPO-NAME\_CLASS-NAME\_H : Code : #endif //
   VOTCA\_VOTCA-REPO-NAME\_CLASS-NAME\_H

-  Never use the "using namespace" in a header file.
-  Avoid using includes in header files. If possible forward declare a
   class instead.

Auto
~~~~

-  Avoid using auto unless the type is very long, the reason being auto
   obscures the underlying type and can make it difficult to discern
   what a variable is meant to be used for.

Classes
~~~~~~~

-  Normally class names in upper case.
-  Order of access modifiers in class definitions should be as follows:
   -  first ``public`` all functions
   -  then ``private``/``protected`` all member variables
   -  then ``private``/``protected`` member functions
-  There is no rule as to where to define a ``public typedef`` in the class.
-  All member variables are ``private``/``protected``.
-  The body of class methods should be placed in a source file or inlined at the end of the header if it exceeds a single line.

Naming in Classes
~~~~~~~~~~~~~~~~~

-  All member variables should be in lower case and end with ``_``.
-  All functions should start with upper case, no ``_`` should exist in their names.
-  Only ``get``/``set`` methods can begin with lower case letters. 
-  For consistency all Ids should start at 0 not 1.

get/set Functions
~~~~~~~~~~~~~~~~~

-  ``get``/``set`` functions should start with a lowercase ``get``/``set`` (these are the only
   functions which should directly ``set``/``get`` a private member variable)
-  ``get`` must return a constant reference and keep the ``class const``:
   ``const int &getId() const;``
-  ``set`` only sets the member, e.g.
   ``void setId(const int &id) { id_ = id; }``

Functions
~~~~~~~~~

-  Functions should remain short.
-  Functions should not have more than one use, so use boolean arguments
   sparingly.

Pointers
~~~~~~~~

-  In general, use pointers sparingly. Most objects are small and a
   copy does not change performance. Use references if you want to avoid copies.
-  If your pointer owns an object (i.e. it has to delete it later) use a
   ``unique_ptr`` to it, so you do not have to call ``delete`` on it
   yourself.
-  If multiple objects own an object and the last object alive should
   delete it, use a ``shared_ptr``.
-  If your object does not have ownership but just wants to visit, you
   can use a raw pointer, but if you can a reference is better.
-  If you ever have to explicitly call ``delete``, you did something
   very wrong.

General
~~~~~~~

-  Do not comment out code, if you do not use it delete it.
-  Variables should have clear and explicit names.
-  Do not duplicate code.
-  Functions should have no more than 3 arguments. Otherwise create a
   class.
-  XYZ positions should be ``Eigen::Vector3d`` from the eigen library.
-  Readability is more important than elegant design.
-  Leave the code better than you found it.
-  Use pointers sparingly and especially try not to pass them around
   objects. Prefer references.
-  Do not write code, which you may use in the future. Only write code
   you will use now. Write code, you need later, later. This avoids
   cluttering the codebase with unused "at some point we will need this
   functions".

VOTCA specifics (indexing, ids, units)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This can all be found in the :doc:`VOTCA_LANGUAGE_GUIDE`.

Testing
-------

Unit Testing
~~~~~~~~~~~~

Each module (tools, csg, xtp) contains a ``src/tests`` folder with the unit
tests of its library.

For every new object and algorithm created there should exist a test. We
use the Boost libraries testing framework, see the
`Boost.Test documentation <https://www.boost.org/doc/libs/release/libs/test/doc/html/index.html>`__.

As an example, take the histogram class ``votca::tools::Histogram``
(``tools/include/votca/tools/histogram.h``):

1. Create a test file ``test_histogram.cc`` in
   `tools/src/tests/ <https://github.com/votca/votca/tree/master/tools/src/tests>`__
   and add its name (without ``.cc``) to the alphabetical list of tests in
   the ``foreach`` of ``tools/src/tests/CMakeLists.txt``. The file has the
   following structure:

   .. code-block:: cpp

       #define BOOST_TEST_MAIN

       #define BOOST_TEST_MODULE histogram_test
       #include <boost/test/unit_test.hpp>

       #include "votca/tools/histogram.h"

       using namespace votca::tools;

       BOOST_AUTO_TEST_SUITE(histogram_test)

       BOOST_AUTO_TEST_CASE(test1) {
         Histogram hist;
         BOOST_CHECK_EQUAL(...);
       }

       BOOST_AUTO_TEST_SUITE_END()

   Replace the '...' with the appropriate checks. For more info on which
   boost test macros to use refer to the boost documentation.

2. Build and run the tests:

   .. code-block:: bash

       cmake --build builddir
       ctest --test-dir builddir -R histogram

CPP Coding Style Guide
-----------------------

VOTCA uses a few auto formatting tools to help enforce the rules.

`clang-format <https://clang.llvm.org/docs/ClangFormat.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Automatically ensures consistent formatting for .cc and .h files. The
style follows the google style formatting rules. Have a look at the
``.clang-format file`` in the `main votca
repository <https://github.com/votca/votca/blob/master/.clang-format>`__
for details.

To run the clang-format function on file.cc.

.. code-block:: bash

    clang-format -i -style=file file.cc

'-i' ensures it will make changes to file.cc, omitting the '-i' will
display the changes without implementing them. '-style=file' ensures the
format is read from the .clang-format file otherwise it will use a
default style guide.

By default tabs should not be used to indent, avoid inserting '\\t', it
is preferable that spaces be used instead.

`autopep8 <https://pypi.org/project/autopep8/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Automatically formats python .py files. We are using the default format
rules of autopep8. To run on file.py and update the file run:

.. code-block:: bash

    autopep8 -i file.py

Automating Formatting
~~~~~~~~~~~~~~~~~~~~~

The above formatters can be automated at every commit using the script
found in the `dev-tools <https://github.com/votca/dev-tools>`__
repository. To use it copy the file ``pre-commit`` to your local .git
subfolder to the hooks folder. E.g.

.. code-block:: bash

    chmod 777 dev-tools/pre-commit  
    cp dev-tools/pre-commit votca/.git/hooks/

The above will make the script executable and then copy it to the local
.git/hooks directory in the votca repository. The script not only
updates the file format of every file staged during a commit it will
also update the license date.

CPP Comment Guide
-----------------

It is preferential that the following guidelines be followed when adding
comments to code:

1. The ``/* */`` comment blocks should be avoided and the ``//`` used in
   their place. This is so that the ``/* */`` comment blocks can be
   easily used for debugging.
2. It would be preferential that the following doxygen commenting
   stencil be used in the header files above each class and function
   description.

   .. code-block:: cpp

       /**
       * \brief function/class summary
       *
       * Detailed function/class description if needed
       *
       * @param[in] - description of parameter 1
       * @param[out] - description of parameter 2
       * @param[in,out] - description of parameter 3
       * :
       * @return - description of return type
       */

Doxygen commenting will help future developers maintain the code, in
its fully compiled state. It may be found at: http://doc.votca.org.

NOTE: Compilation of the doxygen documentation is automated when code is
merged into the :code:`master` votca branch!

Failed Release Builds
---------------------

To prepare votca for distribution on different linux flavors there are
different requirements from the package managers. Some of the
architectures that the package managers support can be quite varied. In
the case that a failure occurs on an architecture, that is not available
to you, there are different approaches for debugging the problem. As an
example, Fedora supports the **ppc64le** architecture. Assuming you have
access to Fedora you can simulate the build process on **ppc64le** with
``mock`` and the Fedora package sources
(`src.fedoraproject.org/rpms/votca <https://src.fedoraproject.org/rpms/votca>`__):

.. code-block:: bash

    dnf install qemu-user-static mock fedpkg
    usermod -a -G mock <username>
    fedpkg clone -a votca && cd votca
    fedpkg srpm
    mock -r fedora-rawhide-ppc64le --forcearch ppc64le votca-*.src.rpm

If the build crashes one can interactively intervene with

.. code-block:: bash

    mock -r fedora-rawhide-ppc64le --forcearch ppc64le --shell

You will also need to install a text editor if you want to change the
source files before running the interactive instance.

.. code-block:: bash

    mock -r fedora-rawhide-ppc64le --forcearch ppc64le --install vim

The same procedure works for other architectures (e.g. **aarch64**) and
other releases (e.g. ``fedora-<version>-<arch>`` or
``epel-<version>-<arch>``; ``mock --list-chroots`` shows them). If you
simply want to build natively, drop ``--forcearch``.
