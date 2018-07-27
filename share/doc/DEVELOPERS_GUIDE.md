
# Developer and Contributor Guide

 The page is designed to give new developers general guidelines for
 implementing code consistent with the VOTCA and cpp style and standard.
 
 - [Reporting Bugs](#reporting-bugs)
 - [CPP Resoures](#cpp-resources)
 - [CPP Tips](#cpp-tips)
 - [Testing](#testing)
 - [CPP Codeing Style Guide](#cpp-codeing-style-guide)
 - [CPP Comment Guide](#cpp-comment-guide)

## Reporting Bugs

 To report a bug please create an issue on the appropriate github repo. Please
 be sure to provide as much information as possible such as:

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
 - [xtp-tutorials](https://github.com/votca/xtp-tutorials/issues)
 - [votca](https://github.com/votca/votca/issues)

## CPP Resources

 A good starting point is to take a look at the cpp standard. Though the code
 has not always consistently followed the cpp standard we now make an
 effort to really enforce it and follow best practices.

 - [Best Practices1](https://www.gitbook.com/book/lefticus/cpp-best-practices/details)
 - [Best Practices2](https://google.github.io/styleguide/cppguide.html)

## CPP Tips

 Here are a few general tips that should be followed:
 
### Files ###
  * each class goes into a separate file
  * filename is name of class in lowercase

### Header Files ###
 * One class, one header.
 * When creating header gruards use the following form, where "VOTCA-REPO-NAME" is replaced by whichever repo the header is in tools/csg/ctp/xtp, and where "CLASS-NAME" is replaced by the name of the class described in the header file:
```
#ifndef _VOTCA_VOTCA-REPO-NAME_CLASS-NAME_H
#define _VOTCA_VOTCA-REPO-NAME_CLASS-NAME_H
:
Code
:
#endif // _VOTCA_VOTCA-REPO-NAME_CLASS-NAME_H
```
 * Never use the "using namespace" in a header file.
 * Avoid using includes in header files. If possible forward declare a class instead.
 
### Braces ###
  * in functions/classes, the { is in the next line
  * for for loops, if, ..., the { is n the same line as if,for

### Classes ###
  * normally begin in upper case
  * exceptions: classes which define "small" types (e.g. `vec`, `matrix`)
  * order in class definition:
    * first `public` all functions
    * then `private`/`protected` all member variables
    * then `private`/`protected` member functions
    * no rule where to define a `public typedef` in the class
  * all member variables are `private`/`public` (again exception e.g. `vec`, `matrix`)
  * maximum one-line-function implementation in class declaration, everything else moves to separate file or inline at end of header.

### Naming in Classes ###
  * all member variables are in lower case and end with `_`
  * all functions start with upper case, no `_` in names
  * exception: `get`/`set` functions
  * for consistency all Ids should start at 0 not 1 

### get/set Functions ###
  * get/set functions start with a lowercase set/get (these are only functions which directly set/get a private member variable)
  * get must return a constant reference and keep the `class const`: `const int &getId() const;`
  * set only sets the member, e.g. `void setId(const int &id) {_id=id;}`
 
### Functions ###
 * Make functions short.
 * Functions should not have more than one use. So use boolean arguments sparingly.
 
### General ###
 * Do not comment out code, if you do not use it delete it.
 * Variables should have clear and explicit names.
 * Do not duplicate code.
 * Functions should have no more than 3 arguments. Otherwise create a class.
 * Do not use raw pointers, `new` and `delete` should not be used, instead use `std::shared_ptr` or `std::unique_ptr`
 * XYZ positions should be described using tools::vec, 3x3 matrices tools::matrix, or classes and or functions in the eigen library.
 * Readability is more important the elegant design.
 * Leave the code better than you found it.
 * Use pointers sparingly and especially try not to pass them around objects. Prefer references.
 * Do not write code, which you may use in the future. Only write code you will use now. Write code, you need later, later. This avoids cluttering the codebase with unused "at some point we will need this functions".
  
## Testing

 Each repository contains a src folder. Within the src folder exists a 
 library folder: libtools, libcsg etc... and a tools folder. A tests folder
 should also exist in the src folder. If it does not you should create one.

 For every new object and algorithm created there should exist a test. We
 use the Boost libraries testing framework. Good documentation can be found 
 here:

 - [Boost link](https://www.ibm.com/developerworks/aix/library/au-ctools1_boost/)
 
 We will outline the general workflow here using the vec object in 
 votca::tools. This object only has a header file it is in:
 tools/include/votca/tools/vec.h
 
 Determine if a tests folder has already been created or not in /src if it
 has not take a look at what was done in the votca-tools repo. 

 1. Create a test file in [tools/src/tests/](https://github.com/votca/tools/tree/master/src/tests)test_vec.cc 
    must have the same name as what appears in the foreach in the 
    CMakeLists.txt file. And place the following contents

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
 
Replace the '...' and ':' with the appropriate syntax. 
For more info on which boost test macros to use refer to the boost documentation

 2. To compile and test the code create a folder tools/build and run the following
    commands: 

```
cmake -DENABLE_TESTING=ON ../
make 
make test
```

Ensure you have an up to date version of cmake or use cmake3 

## CPP Codeing Style Guide

 VOTCA uses the clang formatter to automatically ensure consistent style.
 The style follows the google style fomatting rules. 
 Have a look at `.clang-format file` in the
 [main votca repository](https://github.com/votca/votca/blob/master/.clang-format) for details.


To run the clang-format function on file.cc  

```
clang-format -i -style=file file.cc
```

'-i' ensures it will make change to file.cc, omitting the '-i' will display the
     changes without implementing them.
'-style=file' ensures the format is read from the .clang-format file otherwise 
     it will use a default style guide. 

By default tabs should not be used to indent, avoid inserting '\t', it is
preferable that spaces be used instead. 

## CPP Comment Guide
 
 It is preferential that the following guidelines be followed when adding 
 comments to code:

 1. The `/* */` comment blocks should be avoided and the `//` used in their 
    place. This is so that the `/* */` comment blocks can be easily used for 
    debugging.
 2. It would be preferential that the following doxygen commenting stencil be 
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
