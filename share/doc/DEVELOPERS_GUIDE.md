# Developer and Contributor Guide

The page is designed to give new developers general guidelines for implementing
code consistent with the VOTCA and cpp style and standard.

- [Reporting Bugs](#reporting-bugs)
- [CPP Resoures](#cpp-resources)
- [CPP Tips](#cpp-tips)
- [Testing](#testing)
- [Failed Travis Builds](#failed-travis-builds)
- [CPP Codeing Style Guide](#cpp-codeing-style-guide)
- [CPP Comment Guide](#cpp-comment-guide)
- [Updating Git Submodules](#updating-git-submodules)
- [Merging With Stable](#merging-with-stable)
- [Failed Release Builds](#failed-release-builds)
- [Setting Up A GitLab Runner Server](#gitlab-server)

## Reporting Bugs

To report a bug please create an issue on the appropriate github repo. Please be
sure to provide as much information as possible such as:

- The error messages
- The operating system
- What compiler was used
- What dependencies were installed
- The calculation that was being run

Issues can be directed created on the appropriate github repo:

- [tools](https://github.com/votca/tools/issues)
- [csg](https://github.com/votca/csg/issues)
- [csgapps](https://github.com/votca/csgapps/issues)
- [csg-manual](https://github.com/votca/csg-manual/issues)
- [csg-tutorials](https://github.com/votca/csg-tutorials/issues)
- [ctp](https://github.com/votca/ctp/issues)
- [xtp](https://github.com/votca/xtp/issues)
- [votca](https://github.com/votca/votca/issues)

## CPP Resources

A good starting point is to take a look at the cpp standard. Though the code has
not always consistently followed the cpp standard we now make an effort to
really enforce it and follow best practices.

- [Best Practices1](https://www.gitbook.com/book/lefticus/cpp-best-practices/details)
- [Best Practices2](https://google.github.io/styleguide/cppguide.html)

## CPP Tips

Here are a few general tips that should be followed:

### Files

- each class goes into a separate file
- filename is name of class in lowercase

### Includes

- When including a header file from within the same repo that you are working
  use the relative includes. This consists of using quotation marks i.e.

```
#include "molecule.h"
```

- When including from another repository, for instance you are working in the
  csg repostory and want to include a file from the tools repo use the anglular
  brackets i.e.

```
#include <votca/tools/table.h>
```

### Header Files

- One class, one header.
- When creating header gruards use the following form, where "VOTCA-REPO-NAME"
  is replaced by whichever repo the header is in tools/csg/ctp/xtp, and where
  "CLASS-NAME" is replaced by the name of the class described in the header
  file:

```
#ifndef VOTCA_VOTCA-REPO-NAME_CLASS-NAME_H
#define VOTCA_VOTCA-REPO-NAME_CLASS-NAME_H
:
Code
:
#endif // VOTCA_VOTCA-REPO-NAME_CLASS-NAME_H
```

- Never use the "using namespace" in a header file.
- Avoid using includes in header files. If possible forward declare a class
  instead.

### Braces

- in functions/classes, the { is in the next line
- for for loops, if, ..., the { is n the same line as if,for

### Auto

- avoid using auto unless the type is very long, the reason being auto obscures
  the underlying type and can make it difficult to discern what a variable is
  meant to be used for

### Classes

- normally begin in upper case
- exceptions: classes which define "small" types (e.g. `vec`, `matrix`)
- order in class definition:
  - first `public` all functions
  - then `private`/`protected` all member variables
  - then `private`/`protected` member functions
  - no rule where to define a `public typedef` in the class
- all member variables are `private`/`public` (again exception e.g. `vec`,
  `matrix`)
- maximum one-line-function implementation in class declaration, everything else
  moves to separate file or inline at end of header.

### Naming in Classes

- all member variables are in lower case and end with `_`
- all functions start with upper case, no `_` in names
- exception: `get`/`set` functions
- for consistency all Ids should start at 0 not 1

### get/set Functions

- get/set functions start with a lowercase set/get (these are only functions
  which directly set/get a private member variable)
- get must return a constant reference and keep the `class const`:
  `const int &getId() const;`
- set only sets the member, e.g. `void setId(const int &id) { _id = id; }`

### Functions

- Make functions short.
- Functions should not have more than one use. So use boolean arguments
  sparingly.

### Pointers

- In general, use pointers sparringly. Most objects are small and a copy does
  not change performance. Use references as well
- If your pointer owns an object (i.e. it has to delete it later) use a
  `unique_ptr` to it, so you do not have to call `delete` on it yourself
- If multiple objects own an object and the last object alive should delete it,
  use a `shared_ptr`
- If your object does not have ownership but just wants to visit, you can use a
  raw pointer, but if you can a reference is better.
- If you ever have to explicitly call `delete`, you did something very wrong.

### General

- Do not comment out code, if you do not use it delete it.
- Variables should have clear and explicit names.
- Do not duplicate code.
- Functions should have no more than 3 arguments. Otherwise create a class.
- XYZ positions should be described using tools::vec, 3x3 matrices
  tools::matrix, or classes and or functions in the eigen library.
- Readability is more important the elegant design.
- Leave the code better than you found it.
- Use pointers sparingly and especially try not to pass them around objects.
  Prefer references.
- Do not write code, which you may use in the future. Only write code you will
  use now. Write code, you need later, later. This avoids cluttering the
  codebase with unused "at some point we will need this functions".

## Testing

### Unit Testing

Each repository contains a src folder. Within the src folder exists a library
folder: libtools, libcsg etc... and a tools folder. A tests folder should also
exist in the src folder. If it does not you should create one.

For every new object and algorithm created there should exist a test. We use the
Boost libraries testing framework. Good documentation can be found here:

- [Boost link](https://www.ibm.com/developerworks/aix/library/au-ctools1_boost/)

We will outline the general workflow here using the vec object in votca::tools.
This object only has a header file it is in: tools/include/votca/tools/vec.h

Determine if a tests folder has already been created or not in /src if it has
not take a look at what was done in the votca-tools repo.

1.  Create a test file in
    [tools/src/tests/](https://github.com/votca/tools/tree/master/src/tests)test_vec.cc
    must have the same name as what appears in the foreach in the CMakeLists.txt
    file. And place the following contents

```
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
```

Replace the '...' and ':' with the appropriate syntax. For more info on which
boost test macros to use refer to the boost documentation

2.  To compile and test the code create a folder tools/build and run the
    following commands:

```
cmake -DENABLE_TESTING=ON ../
make
make test
```

Ensure you have an up to date version of cmake or use cmake3

### Testing Across Repos

There may come a case where changes have to be committed across more than one
repo at the same time. Attempting to merge one repo at a time will cause the
continuous integration to fail as changes in the other repos will not be pulled
in. To do this correctly the following steps should be taken.

Assuming you are in the votca/votca repository:

    git checkout <base_branch>
    git submodule update
    git checkout -b <some_descriptive_branch_name>
    git submodule foreach git remote update
    git -C <module1> checkout <sha_or_branch_of_module1_to_test>
    git -C <module2> checkout <sha_or_branch_of_module2_to_test>
    git add <module1> <module2>
    git commit -m "test <module1> with <module2>"
    git push origin <some_descriptive_branch_name>

1. Here `base_branch` will typically be the master or stable branch.
```
git checkout <base_branch>
```
2. The submodules are updated to be sure they have incorporated the latest
   changes in your local repository
```
git submodule update
```
3. Create a branch with a descriptive name
```
git checkout -b <some_descriptive_name>
```
4. Update each of the submodules, by pulling in any remote changes to the
   submodules.
```
git submodule foreach git remote update
```
5. '-C' changes directory to the submodule directory and then checks out the
   appropriate commit
```
git -C <module1> checkout <sha_or_branch_of_module1_to_test>  
git -C <module2> checkout <sha_or_branch_of_module2_to_test>
```
6. The changes are then added and commited
```
git add <module1> <module2>  
git commit -m "test <module1> with <module2>"
```
7. Finally, they are pushed to the remote branch
```
git push origin <some_descriptive_branch_name>
```
A pull request is then made for the votca/votca repo using the branch name. Once
the branch passes all tests it can be merged. Pull requests for each of repos
changed can then be made. They will now compile against the updated votca/votca
repo. Once they pass their tests they can be merged. If a pull request was
already made the travis tests may simply need to be restarted.

## Failed Travis Builds

There may come a time where one of the docker builds fails. It may be the case
that the error message is clear and it can be reproduced on your host os.
However, in the case that the error is specific to the enviorment used in the
build the local enviornment can be simulated using a docker container.

Before you can use this approach docker must be installed on your host OS. Begin
by running a docker image the default is:

```
docker run -it votca/buildenv:fedora /bin/bash
```

This will run an interative docker container which you can interact with in bash
. The next commands will need to be adjusted to whatever local environment you
need to reproduce to test the error in the travis build.

## CPP Codeing Style Guide

VOTCA uses the clang formatter to automatically ensure consistent style. The
style follows the google style fomatting rules. Have a look at
`.clang-format file` in the
[main votca repository](https://github.com/votca/votca/blob/master/.clang-format)
for details.

To run the clang-format function on file.cc

```
clang-format -i -style=file file.cc
```

'-i' ensures it will make change to file.cc, omitting the '-i' will display the
changes without implementing them. '-style=file' ensures the format is read from
the .clang-format file otherwise it will use a default style guide.

By default tabs should not be used to indent, avoid inserting '\t', it is
preferable that spaces be used instead.

Clang formatting can be automated at every commit using the script found in the
[dev-tools](https://github.com/votca/dev-tools) repository. To use it copy the
file `pre-commit` to your local .git subfolder to the hooks folder. E.g.   
 chmod 777 dev-tools/pre-commit  
 cp dev-tools/pre-commit tools/.git/hooks/

The above will make the script executable and then copy it to the local
.git/hooks directory in the tools repository. The script not only updates the
file format of every file staged during a commit it will also update the license
date.

## CPP Comment Guide

It is preferential that the following guidelines be followed when adding
comments to code:

1.  The `/* */` comment blocks should be avoided and the `//` used in their
    place. This is so that the `/* */` comment blocks can be easily used for
    debugging.
2.  It would be preferential that the following doxygen commenting stencil be
    used in the header files above each class and function description.

```
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
```

The doxygen commenting will help future developers maintain the code, in its
fully compiled state it may be found at: http://doc.votca.org

NOTE: Compilation of the doxygen documentation is automated when code is merged
into the master votca branch!

## Updating Git Submodules

Votca with all of its repos can be build by using the parent
[votca repo](https://github.com/votca/votca). All the other necessary repos
appear as submodules in the parent repo. It is worth noting that the submodules
are not automatically updated whenever changes are made to their respective
master branches. In essence a submodule refers to a specific commit of the repo
it represents. If a new commit is merged into the master branch of a repository
the submodule state in the parent repo has to be updated for the commit to
propagate to the parent votca repository.

To update the state of a submodule the following commands can be used:

```
git submodule foreach git checkout master
git submodule foreach git pull
git add -u
git commit -m "update all submodules"
```

## Merging With Stable

When creating a pull request to merge a branch with master the ci test build
will by default build agains the master votca branch. This is fine if you are
merging with the master branch of a given repo. However, if you would like to
merge a bug fix to stable it becomes a problem because you want to build against
the stable branch.

Thus to merge with stable the name of the branch to be merged must be prepended
with **for/stable/**. As an example if I have a bug fix for stable I would name
it **for/stable/bug-fix-1-because-I-am-awesome**. When the ci runs it should
proceed to make the build correctly.

## Failed Release Builds

To prepare votca for distribution on different linux flavors there are different
requirements from the package managers. Some of the architectures that the
package managers support can be quite varied. In the case that a failure occurs
on an architecture that is not available to you there are different approaches
to debugging the problem. As an example fedora dnf has extended support to the
**pcc64le** architecture. Assuming you have access to fedora you can run the
following commands to simulate the build process on the **pcc64le**
architecture:

```
dnf update
dnf install qemu-user-static dnf-utils
usermod -a -G mock <username>
mock -r epel-7-ppc64le --forcearch ppc64le --dnf --init
wget https://raw.githubusercontent.com/votca/fedora-copr/master/votca.spec
spectool -g votca.spec
rpmbuild -D"_sourcedir ${PWD}" -D"_srcrpmdir ${PWD}" -bs votca.spec
mock -r epel-7-ppc64le --forcearch ppc64le --dnf --no-clean votca-1.5-1.*.src.rpm
```

Here, votca-1.5-1 should be replaced with the correct version. The above
commands would setup and run the dnf installation process on the **pcc64le**
enviroment. If a bug was found and the build crashes one can interactively
intervene by issuing the following command:

```
mock -r epel-7-ppc64le --forcearch ppc64le --shell
```

You will also need to install a text editor if you want to change the source
files before running the interactive instance.

```
mock -r epel-7-ppc64le --forcearch ppc64le --install vim
```

Note: we have used this process with the **ppc64le** architecture as an example,
but the same procedure can be extended with different architectures and diferent
operating systems. For example you could use the **aarch64** or **armv7hl**
architecture in place of **pcc64le**. You could also replace the
**epel-7-ppc64le** os-architecure to **fedora-28-ppc64le**,
**fedora-27-aarch64** or some other combination. A final point, if you simply
want to build natively for instance if you are running fedora on an **x86_64**
machine the `frocearch pcc64le` in the above case could just be dropped.

## Gitlab Server

To setup a gitlab server the instructions on the official
[gitlab website](https://docs.gitlab.com/runner/register/) should be followed.
To obtain the access token contact one of the votca administrators. It is also
advisable to setup a cron job to help monitor the docker containers, images and
volumes which can use a substantial amount of space if not cleaned. A docker
cleanup script has been added in the
[dev-tools](https://github.com/votca/dev-tools) repo, instructions on how to set
it up are provided in the script.
