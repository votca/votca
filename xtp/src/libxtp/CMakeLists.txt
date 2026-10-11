file(GLOB_RECURSE VOTCA_SOURCES *.cc *.cpp)

get_filename_component(_PARENT_PATH "${CMAKE_CURRENT_BINARY_DIR}" DIRECTORY)

# The translation units that include libint2.hpp need several GB each to
# compile (2-3.4 GB with libint2 2.7, 5-7 GB with the libint2 of Fedora
# rawhide in 2026), almost all of it for parsing the header. Several of them
# at once exhaust the memory of a 16 GB machine under a parallel build. They
# are built as a separate object library in a job pool of one, so they never
# run concurrently; the other files keep the full parallelism. The pool is
# also used by the unit test that includes libint2.hpp. JOB_POOL_COMPILE is
# honoured by Ninja and ignored by the other generators.
set(VOTCA_LIBINT2_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/libint2_calls.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/libint2_derivative_calls.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${VOTCA_LIBINT2_SOURCES})
set_property(GLOBAL APPEND PROPERTY JOB_POOLS votca_libint2=1)

add_library(votca_xtp_libint2 OBJECT ${VOTCA_LIBINT2_SOURCES})
target_compile_features(votca_xtp_libint2 PRIVATE cxx_std_17)
set_target_properties(votca_xtp_libint2 PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  JOB_POOL_COMPILE votca_libint2)

add_library(votca_xtp  ${VOTCA_SOURCES} $<TARGET_OBJECTS:votca_xtp_libint2>)
target_compile_features(votca_xtp PUBLIC cxx_std_17)

set_target_properties(votca_xtp PROPERTIES SOVERSION ${SOVERSION})
set_target_properties(votca_xtp PROPERTIES VOTCA_SHARE ${PROJECT_SOURCE_DIR}/xtp/share)
set_target_properties(votca_xtp PROPERTIES XTP_PATH ${_PARENT_PATH}/tools)

# to treat libint as system to not trigger -Werror for library internal warnings
add_library(Libint2_system INTERFACE)
target_link_libraries(Libint2_system INTERFACE Libint2::int2)

get_target_property(_libint_incs Libint2::int2 INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(Libint2_system SYSTEM INTERFACE ${_libint_incs})

target_link_libraries(votca_xtp PUBLIC VOTCA::votca_csg VOTCA::votca_tools Boost::boost Eigen3::Eigen ${HDF5_LIBRARIES} PRIVATE Libxc::xc Libint2_system ECPINT::ecpint Boost::program_options Boost::timer openbabel)
# same usage requirements for the libint2 objects (an object library only
# takes over include directories, definitions and options from these)
target_link_libraries(votca_xtp_libint2 PRIVATE VOTCA::votca_csg VOTCA::votca_tools Boost::boost Eigen3::Eigen Libxc::xc Libint2_system ECPINT::ecpint Boost::program_options Boost::timer openbabel)

# OpenBabel3's own "openbabel" target is an older-style, unnamespaced
# target (confirmed directly, from OpenBabel's own official CMake
# example, doc/examples/static_executable/CMakeLists.txt, which links
# against it this same way) -- unlike a modern imported target, it may
# not automatically carry its own include directories, so these are
# added explicitly here too, matching that same official example's own
# separate include_directories(${OpenBabel3_INCLUDE_DIRS}) call.
#
# Real, direct, genuine fix, confirmed directly from a real CI compile
# failure (Ubuntu, gcc 16): SYSTEM added here -- the exact same,
# already-established pattern this same file already uses, right
# above, for Libint2_system, and for the exact same real reason ("to
# treat libint as system to not trigger -Werror for library internal
# warnings") -- OpenBabel3's own real headers (atom.h/base.h/
# forcefield.h/plugin.h) trigger multiple, real -Werror-promoted
# warnings of their own (unused-parameter/conversion/
# deprecated-declarations/shadow) that are entirely OpenBabel3's own,
# not this project's own code at all, and genuinely out of this
# project's own control to fix directly.
target_include_directories(votca_xtp SYSTEM PRIVATE ${OpenBabel3_INCLUDE_DIRS})
target_include_directories(votca_xtp_libint2 SYSTEM PRIVATE ${OpenBabel3_INCLUDE_DIRS})

# PROJECT_BINARY_DIR/xtp/include/votca/xtp for votca_xtp_config.h
target_include_directories(votca_xtp PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/xtp/include> ${HDF5_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/xtp/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/xtp/include/votca/xtp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(votca_xtp_libint2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/xtp/include ${HDF5_INCLUDE_DIRS}
  ${PROJECT_BINARY_DIR}/xtp/include
  ${PROJECT_BINARY_DIR}/xtp/include/votca/xtp)
install(TARGETS votca_xtp EXPORT VOTCA_XTP_Targets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT VOTCA_XTP_Targets FILE VOTCA_XTP_Targets.cmake NAMESPACE VOTCA:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VOTCA_XTP)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("VOTCA_XTPConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
configure_file(VOTCA_XTPConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_XTPConfig.cmake" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_XTPConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_XTPConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VOTCA_XTP)
add_library(VOTCA::votca_xtp ALIAS votca_xtp)
