########################################################################
#Find external packages
########################################################################
find_program(GMX_EXECUTABLE NAMES gmx_d gmx)
find_package_handle_standard_args(GMX REQUIRED_VARS GMX_EXECUTABLE)

if(BUILD_TESTING)
  set(INTEGRATIONTEST_TOLERANCE "5e-5" CACHE STRING "Tolerance for the integration tests")
  mark_as_advanced(INTEGRATIONTEST_TOLERANCE)
  find_package(UnixCommands REQUIRED)
  set_package_properties(UnixCommands PROPERTIES TYPE OPTIONAL PURPOSE "Used in testing")
endif(BUILD_TESTING)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFileCXX)
foreach(HEADER algorithm cassert fstream iomanip iostream list map numeric sstream stdexcept string vector cstdlib)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

######################################
# Include the following subdirectory # 
######################################
add_subdirectory(src)
add_subdirectory(include/votca/csg)
add_subdirectory(share)
