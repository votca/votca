For more detailed information about the changes see the history of the [repository](https://github.com/votca/tools/commits/stable).

## Version 1.6-dev
* replaced fftw with Eigen3 fft 
* replaced random with std::random implementation

## Version 1.5.1 (released 20.11.19)
 * switch votca_compare to python3
 * more clear error message for xmlreader
 * fix build on CentOs7

## Version 1.5 _SuperVictor_ (released 31.01.19)

* enable gitlab CI

## Version 1.5_rc3 (released 19.01.19)

* travis: fixed bug in building tags

## Version 1.5_rc2 (released 16.01.19)

* histogramnew: refactor to fix memory lead

## Version 1.5_rc1 (released 28.12.18)

* added unit testing 
* Eigen added as a dependency
* removed GSL dependency

## Version 1.2.2

* cmake: allow static fftw and gsl
* added dummy c function for cmake
* fixed conflicting type headers (real was defined)

## Version 1.2.1

* fixed soname of libs
* improved cmake checks and error messages
* fixed pkg-config file

## Version 1.2 (SuperDoris)

* changed buildsystem to cmake
* added database class through sqlite3

## Version 1.1.2

* fixed head of votca-tools manpage

## Version 1.1.1

* just a bugfix release

## Version 1.1 (SuperAnn)

* added thread class
* added spline class, with akima spline, linear spline
* added man7 page
* random.cc: avoid calling of exit()
* added lexical cast class

## Version 1.0.1

* fixed a bug in VOTCARC.csh for empty LD_LIBRARY_PATH
* completion file has moved back to csg
* added --disable-rc-files to configure
* updated bundled libtool to 2.2.10

## Version 1.0

* added libvotca_expat to allow compiling without expat
* allow comments in tables
* added application class to create standardized applications
* all boost dependecy are now in tools
* fixes in table format, flags is always last row now
* allow compling without fftw (needed for csg_boltzmann only)
* allow compling without gsl (needed for csg_resample and csg_fmatch)

## Version 1.0_rc5

* using libexpat instead of libxml2
* added libvotca_boost to allow compiling without boost
* using pkg-config to detect package flags
* compiles under AIX with xlC
* added VOTCARC to initialize all votca parts
* updated configure, see --help

## Version 1.0_rc4

* using libtool to build shared libs

## Version 1.0_rc3

* corrected bug in tokenizer
* fixed a bug in calculation of version string
* some fixes concerning autotools

## Version 1.0_rc2

* fixed typo in calculation of version string
* added NOTICE and LICENSE to the dist tarball

## Version 1.0_rc1

* initial version
