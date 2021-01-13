Developer and Contributor Guide
===============================

This page is designed to give new developers general guidelines for
implementing code consistent with the VOTCA and cpp style and standard.

-  `Reporting Bugs <#reporting-bugs>`__
-  `Making a release <#making-a-release>`__
-  `CPP Resources <#cpp-resources>`__
-  `CPP Coding Rules <#CPP-Coding-Rules>`__
-  `Testing <#testing>`__
-  `CPP Coding Style Guide <#cpp-coding-style-guide>`__
-  `CPP Comment Guide <#cpp-comment-guide>`__
-  `Updating Git Submodules <#updating-git-submodules>`__
-  `Updates from stable <#updates-from-stable>`__
-  `Failed Release Builds <#failed-release-builds>`__

Reporting Bugs
--------------

To report a bug, please create an issue on the appropriate github repo.
Please be sure to provide as much information as possible such as:

-  The error messages
-  The operating system
-  What compiler was used
-  What dependencies were installed
-  The calculation that was being run

Issues can be directly created on the appropriate github repo:

-  `tools <https://github.com/votca/tools/issues>`__
-  `csg <https://github.com/votca/csg/issues>`__
-  `csg-tutorials <https://github.com/votca/csg-tutorials/issues>`__
-  `xtp <https://github.com/votca/xtp/issues>`__
-  `votca <https://github.com/votca/votca/issues>`__

Formatting code
---------------

VOTCA uses ``clang-format`` to format the code, code that is not
properly formatted is automatically rejected. The style files can be
found in each repo. If you have a version ``>7.1`` of ``clang-format``
you can just run ``make format`` and your code will be formatted. For an even
more automated version, see `VOTCA dev-tools <#votca-dev-tools>`__

VOTCA dev-tools
---------------

Running clang-format on every commit can be a drag, as can changing the
copyright in every header. Building artifacts locally from a Gitlab run
also takes multiple steps. Fortunately, you will find small scripts in the
`dev-tools repo <https://github.com/votca/dev-tools>`__, which can
automate this.

VOTCA Continuous Integration (Github Actions)
---------------------------------------------

Each pull request to master in the tools, csg, csg-tutorials, xtp, xtp-tutorials or votca repository 
is built on a machine in the cloud using `Github actions <https://docs.github.com/en/actions>`__ (There is still some Gitlab for the GPU builds).

VOTCA can be built on various linux distributions, which are not all natively supported by Github actions. For non natively supported distributions, 
instead of using the default virtual machines, VOTCA first builds and then runs a `docker container <https://www.docker.com/resources/what-container>`__ for each Pull Request. The container contains all the necessary dependencies of VOTCA (see :code:`buildenv` below)

The docker images can be found at `Docker Hub <https://hub.docker.com/u/votca>`__. The **votca/buildenv** containers are the basic containers, which contain all the dependencies VOTCA requires; VOTCA code itself is not included. The **votca/buildenv** can be found on `VOTCA's GitHub Container registry <https://github.com/orgs/votca/packages>`__. 
The actual containers used for running the test builds are built on top of the **votca/buildenv** containers, the resulting **votca/votca** container can be found on `Docker Hub <https://hub.docker.com/u/votca>`__ as well as `VOTCA's GitHub Container registry <https://github.com/orgs/votca/packages>`__.

More information can be found in the `Github workflow files <https://github.com/votca/votca/tree/master/.github/workflows>`__.

Making a Release
----------------

Similar to the VOTCA containers, releases are also handled by Github actions. :code:`votca/votca` has a :code:`release` workflow that can only be triggered manually.
To trigger it go `here <https://github.com/votca/votca/actions?query=workflow%3Arelease>`_. The release can only be made from the 
:code:`stable` branch, but testing the creation of a release can be triggered on any branch. To make a release, trigger the action from the
:code:`stable` branch, pick a new release tag in the :code:`release tag` box (all CHANGELOG files should already contain a section with the tag, but the date will be updated) and type :code:`yesyesyes` into the deploy box. A new release will trigger the creation of the release tag in all involved submodules (plus pull requests for the :code:`stable` to :code:`master` branch, see `below <#updates-from-stable>`__). 

CPP Resources
-------------

A good starting point, is to take a look at the cpp standard. Though the
code has not always consistently followed the cpp standard we now make
an effort to really enforce it and follow best practices.

-  `Best
   Practices1 <https://www.gitbook.com/book/lefticus/cpp-best-practices/details>`__
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

-  When including from another repository, for instance you are working
   in the csg repostory and want to include a file from the tools repo
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
-  All member variables are ``private``/``public``.
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
   ``void setId(const int &id) { _id = id; }``

Functions
~~~~~~~~~

-  Functions should remain short.
-  Functions should not have more than one use, so use boolean arguments
   sparingly.

Pointers
~~~~~~~~

-  In general, use pointers sparringly. Most objects are small and a
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

This can all be found here `VOTCA\_LANGUAGE\_GUIDE <share/doc/VOTCA_LANGUAGE_GUIDE.rst>`__

Testing
-------

Unit Testing
~~~~~~~~~~~~

Each repository contains a src folder. Within the src folder exists a
library folder: libtools, libcsg etc... and a tools folder. A tests
folder should also exist in the src folder. If it does not you should
create one.

For every new object and algorithm created there should exist a test. We
use the Boost libraries testing framework. Good documentation can be
found here:

-  `Boost
   link <https://www.ibm.com/developerworks/aix/library/au-ctools1_boost/>`__

We will outline the general workflow here using the vec object in
votca::tools. This object only has a header file it is in:
tools/include/votca/tools/vec.h.

Determine if a tests folder has already been created or not in /src. If
it has not, take a look at what was done in the votca-tools repo.

1. Create a test file in
   `tools/src/tests/ <https://github.com/votca/tools/tree/master/src/tests>`__\ test\_vec.cc
   must have the same name as what appears in the foreach in the
   CMakeLists.txt file. And place the following contents:

   ::

       #define BOOST_TEST_MAIN

       #define BOOST_TEST_MODULE vec_test
       #include <boost/test/unit_test.hpp>
       #include <exception>

       #include <votca/tools/vec.h>

       using namespace std;
       using namespace votca::tools;

       BOOST_AUTO_TEST_SUITE(vec_test)


       BOOST_AUTO_TEST_CASE(test1){
         vecv;
         BOOST_CHECK_EQUAL(...);
         BOOST_CHECK_EQUAL(...);
         :
       }
       BOOST_AUTO_TEST_CASE(test2){
         vecv;
         BOOST_CHECK_EQUAL(...);
         BOOST_CHECK_EQUAL(...);
         :
       }
       :
       BOOST_AUTO_TEST_SUITE_END()

Replace the '...' and ':' with the appropriate syntax. For more info on
which boost test macros to use refer to the boost documentation

2. To compile and test the code create a folder tools/build and run the
   following commands:

   ::

       cmake -DENABLE_TESTING=ON ../
       make
       make test

Ensure you have an up to date version of cmake or use cmake3.

Testing Across Repos
~~~~~~~~~~~~~~~~~~~~

There may come a time where changes have to be committed across more
than one repo at the same time. Attempting to merge one repo at a time
will cause the continuous integration to fail as changes in the other
repos will not be pulled in. To do this correctly the following steps
should be taken.

Assuming you are in the votca/votca repository:

::

    git checkout <base_branch>
    git submodule update
    git checkout -b <some_descriptive_branch_name>
    git submodule foreach git remote update
    git -C <module1> checkout <sha_or_branch_of_module1_to_test>
    git -C <module2> checkout <sha_or_branch_of_module2_to_test>
    git add <module1> <module2>
    git commit -m "test <module1> with <module2>"
    git push origin <some_descriptive_branch_name>

1. Here ``base_branch`` will typically be the :code:`master` or :code:`stable` branch.

   ::

       git checkout <base_branch>

2. The submodules are updated to be sure they have incorporated the
   latest changes in your local repository.

   ::

       git submodule update

3. Create a branch with a descriptive name.

   ::

       git checkout -b <some_descriptive_name>

4. Update each of the submodules, by pulling in any remote changes to
   the submodules.

   ::

       git submodule foreach git remote update

5. '-C' changes directory to the submodule directory and then checks out
   the appropriate commit.

   ::

       git -C <module1> checkout <sha_or_branch_of_module1_to_test>  
       git -C <module2> checkout <sha_or_branch_of_module2_to_test>

6. The changes are then added and commited.

   ::

       git add <module1> <module2>  
       git commit -m "test <module1> with <module2>"

7. Finally, they are pushed to the remote branch.

   ::

       git push origin <some_descriptive_branch_name>

A pull request is then made for the votca/votca repo using the branch
name. Once the branch passes all tests, it can be merged. Pull requests
for each of the repos changed can then be made. They will now compile
against the updated votca/votca repo. Once they pass their tests, they
can be merged. If a pull request was already made, the travis tests may
simply need to be restarted.

CPP Coding Style Guide
-----------------------

VOTCA uses a few auto formatting tools to help enforce the rules.

`clang-format <https://clang.llvm.org/docs/ClangFormat.html>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Automatically ensures consistent formatting for .cc and .h files. The
style follows the google style fomatting rules. Have a look at the
``.clang-format file`` in the `main votca
repository <https://github.com/votca/votca/blob/master/.clang-format>`__
for details.

To run the clang-format function on file.cc.

::

    clang-format -i -style=file file.cc

'-i' ensures it will make changes to file.cc, omitting the '-i' will
display the changes without implementing them. '-style=file' ensures the
format is read from the .clang-format file otherwise it will use a
default style guide.

By default tabs should not be used to indent, avoid inserting '\\t', it
is preferable that spaces be used instead.

`autopep8 <https://pypi.org/project/autopep8/0.8/>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Automatically formats python .py files. We are useing the default format
rules of autopep8. To run on file.py and update the file run:

::

    autopep8 -i file.py

Automating Formatting
~~~~~~~~~~~~~~~~~~~~~

The above formatters can be automated at every commit using the script
found in the `dev-tools <https://github.com/votca/dev-tools>`__
repository. To use it copy the file ``pre-commit`` to your local .git
subfolder to the hooks folder. E.g.

::

    chmod 777 dev-tools/pre-commit  
    cp dev-tools/pre-commit tools/.git/hooks/

The above will make the script executable and then copy it to the local
.git/hooks directory in the tools repository. The script not only
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

   ::

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

Updating Git Submodules
-----------------------

VOTCA with all of its repos can be build by using the parent `votca
repo <https://github.com/votca/votca>`__. All the other necessary repos
appear as submodules in the parent repo. It is worth noting, the
submodules are automatically updated through a pull request whenever changes are made to
their respective :code:`master` branches. In essence, a submodule refers to a
specific commit of the repo it represents. 

Normally, it is not necessary, but occassionally a new commit must be manually
merged into the :code:`master` branch of a child repository. If this occurs, the
submodule state in the parent repo also has to be updated to reflect the latest
commit of the child repo. 

To update the state of a submodule the following commands can be used:

::

    git submodule foreach git checkout master
    git submodule foreach git pull
    git add -u
    git commit -m "update all submodules"


Updates from :code:`stable`

The :code:`stable` branch contains the latest release with the most uptodate bug fixes since the release.
Only in very limited circumstances should new features be merged into the :code:`stable` branch.
Developers can add bug fixes by making a pull request with the :code:`stable` branch as target.

As the :code:`master` branch of each repo is a child of each repo's :code:`stable` branch,  
any bugfix added to a repos :code:`stable` branch will also need to be pulled into its :code:`master` branch. If the bugfix
is added in one of the child repositories (not :code:`votca/votca`) then :code:`votca/votca` will also need to
reflect these changes.

Keeping the repositories synchronised can be difficult. In order to help keep the :code:`master` branches and :code:`votca/votca`
synchronised with changes in the :code:`stable` branch of a child repository the generation of four pull requests are
automatically generated anytime a bugfix is made to the :code:`stable` branch of a child repository.

E.g. if :code:`hot-bug-fix` is merged into the :code:`stable` branch of :code:`tools`:

1. A pull request is created to merge :code:`stable` from :code:`tools` (child repo) into :code:`master` of :code:`tools` (child repo).
2. A pull request is created to merge :code:`stable` from :code:`tools` (child repo) into :code:`stable` of :code:`votca/votca` (parent repo). This shoud consiste of updating the submodules in the :code:`stable` branch of votca/votca.
3. A pull request is created to merge :code:`master` from :code:`tools` (child repo) into :code:`master` of :code:`votca/votca` (parent repo). Again this should consist of updating the submodules but in the :code:`master` branch of votca/votca.
4. Finally, a pull request is made from :code:`stable` from :code:`votca/votca` (parent repo) to :code:`master` of :code:`votca/votca` (parent repo).

To minimize manual work, it is usually best to merge the pull requests in the order that hey have been shown in the example. 

Failed Release Builds
---------------------

To prepare votca for distribution on different linux flavors there are
different requirements from the package managers. Some of the
architectures that the package managers support can be quite varied. In
the case that a failure occurs on an architecture, that is not available
to you, there are different approaches for debugging the problem. As an
example, fedora dnf has extended support to the **pcc64le** architecture.
Assuming you have access to fedora you can run the following commands to
simulate the build process on the **pcc64le** architecture:

::

    dnf update
    dnf install qemu-user-static dnf-utils
    usermod -a -G mock <username>
    mock -r epel-7-ppc64le --forcearch ppc64le --dnf --init
    wget https://raw.githubusercontent.com/votca/fedora-copr/master/votca.spec
    spectool -g votca.spec
    rpmbuild -D"_sourcedir ${PWD}" -D"_srcrpmdir ${PWD}" -bs votca.spec
    mock -r epel-7-ppc64le --forcearch ppc64le --dnf --no-clean votca-1.5-1.*.src.rpm

Here, votca-1.5-1 should be replaced with the correct version. The above
commands would setup and run the dnf installation process on the
**pcc64le** enviroment. If a bug was found and the build crashes one can
interactively intervene by issuing the following command:

::

    mock -r epel-7-ppc64le --forcearch ppc64le --shell

You will also need to install a text editor if you want to change the
source files before running the interactive instance.

::

    mock -r epel-7-ppc64le --forcearch ppc64le --install vim

Note: we have used this process with the **ppc64le** architecture as an
example, but the same procedure can be extended with different
architectures and diferent operating systems. For example, you could use
the **aarch64** or **armv7hl** architecture in place of **pcc64le**. You
could also replace the **epel-7-ppc64le** os-architecure to
**fedora-28-ppc64le**, **fedora-27-aarch64** or some other combination.
A final point, if you simply want to build natively, for instance if you
are running fedora on an **x86\_64** machine, the ``frocearch pcc64le``
in the above case could just be dropped.
