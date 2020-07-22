For more detailed information about the changes see the history of the
`repository <https://github.com/votca/tools/commits/stable>`__.

Version 1.7-dev
===============

-  add support of txt2tqs 3.4 (#208)
-  fix MKL:: target creation (#207, #265, #269, #270)
-  use local includes (#204)
-  add more support for units (#200)
-  add methods to update user input with defaults (#222)
-  fix usage of PROJECT\_SOURCE\_DIR in cmake (#230)
-  fix merging calculator defaults (#226, #227)
-  document cmake options and packages (#231)
-  make fftw3 a public dependency (#235)
-  rename WITH\_RC\_FILES to INSTALL\_RC\_FILES (#236)
-  check calculator input options (#232, #233)
-  allow calculator choices to be a list (#239)
-  inject defaults into calculator values (#241)
-  fixed bug introduced by defaults injection (#241, #242)
-  make votca\_compare always executable in builddir (#244)
-  add rel. compare to votca\_compare (#143)
-  standardize header formatting (#246)
-  move CI to GitHub Ations (#251, #254, #256, #258, #260,
   #274, #275)
-  improve mkl detection (#257)
-  clean up installed cmake find modules (#263)
-  update codacy badge (#264)
-  add mkl builds to CI (#267, #271)
-  convert markdown files to rst files (#272)
-  enable CXX only in CMake (#273)
-  format code with clang-10 (#277)
-  fix clang-9 conversion warning (#276)

Version 1.6.2 (released XX.07.20)
=================================

-  format code with clang-10 (#278)
-  move CI to GitHub Actions (#280, #285)
-  fix build with mkl (#283)

Version 1.6.1 (released 21.06.20)
=================================

-  fix build with mkl (#229)
-  fix build with non-system libfftw (#234)
-  fix CI on Ubuntu-20.04 (#237)
-  fix bug related to calling c abs instead of c++ fabs (#248)
-  updated floating point comparison in akimaspline.h (#248)
-  fix compile error in structure parameters by adding hash function
   (#248)

Version 1.6 *SuperPelagia* (released 17.04.20)
==============================================

-  fix clang-10 warnings (#217)
-  clean up VOTCARC.\* (#220)
-  fix exported target with mkl (#223, #224)

Version 1.6\_rc2 (released 10.02.20)
====================================

-  fix inheritance bug in cmake files (#197)
-  fix bug in table class (#196)
-  fix build on 32-bit archs (#199)
-  remove old doxygen target (#211)

Version 1.6\_rc1 (released 04.12.19)
====================================

-  replaced fftw with Eigen3 fft
-  replaced random with std::random implementation
-  added many unit tests
-  formated code with clang-format
-  major cmake refactor

Version 1.5.1 (released 20.11.19)
=================================

-  switch votca\_compare to python3
-  more clear error message for xmlreader
-  fix build on CentOs7

Version 1.5 *SuperVictor* (released 31.01.19)
=============================================

-  enable gitlab CI

Version 1.5\_rc3 (released 19.01.19)
====================================

-  travis: fixed bug in building tags

Version 1.5\_rc2 (released 16.01.19)
====================================

-  histogramnew: refactor to fix memory lead

Version 1.5\_rc1 (released 28.12.18)
====================================

-  added unit testing
-  Eigen added as a dependency
-  removed GSL dependency

Version 1.2.2
=============

-  cmake: allow static fftw and gsl
-  added dummy c function for cmake
-  fixed conflicting type headers (real was defined)

Version 1.2.1
=============

-  fixed soname of libs
-  improved cmake checks and error messages
-  fixed pkg-config file

Version 1.2 (SuperDoris)
========================

-  changed buildsystem to cmake
-  added database class through sqlite3

Version 1.1.2
=============

-  fixed head of votca-tools manpage

Version 1.1.1
=============

-  just a bugfix release

Version 1.1 (SuperAnn)
======================

-  added thread class
-  added spline class, with akima spline, linear spline
-  added man7 page
-  random.cc: avoid calling of exit()
-  added lexical cast class

Version 1.0.1
=============

-  fixed a bug in VOTCARC.csh for empty LD\_LIBRARY\_PATH
-  completion file has moved back to csg
-  added --disable-rc-files to configure
-  updated bundled libtool to 2.2.10

Version 1.0
===========

-  added libvotca\_expat to allow compiling without expat
-  allow comments in tables
-  added application class to create standardized applications
-  all boost dependecy are now in tools
-  fixes in table format, flags is always last row now
-  allow compling without fftw (needed for csg\_boltzmann only)
-  allow compling without gsl (needed for csg\_resample and csg\_fmatch)

Version 1.0\_rc5
================

-  using libexpat instead of libxml2
-  added libvotca\_boost to allow compiling without boost
-  using pkg-config to detect package flags
-  compiles under AIX with xlC
-  added VOTCARC to initialize all votca parts
-  updated configure, see --help

Version 1.0\_rc4
================

-  using libtool to build shared libs

Version 1.0\_rc3
================

-  corrected bug in tokenizer
-  fixed a bug in calculation of version string
-  some fixes concerning autotools

Version 1.0\_rc2
================

-  fixed typo in calculation of version string
-  added NOTICE and LICENSE to the dist tarball

Version 1.0\_rc1
================

-  initial version
