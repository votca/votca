cmake_minimum_required(VERSION 3.10)

project(xtp-tutorials LANGUAGES NONE)

set(PROJECT_VERSION "1.7-dev")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)
include(GNUInstallDirs)
include(FeatureSummary)

find_package(VOTCA_XTP REQUIRED NO_MODULE)
set_package_properties(VOTCA_XTP PROPERTIES TYPE REQUIRED PURPOSE "Votca xtp library")
if(DEFINED VOTCA_XTP_VERSION AND NOT VOTCA_XTP_VERSION STREQUAL PROJECT_VERSION)
  message(FATAL_ERROR "Incompatible VOTCA XTP version found (needed ${PROJECT_VERSION}, found ${VOTCA_XTP_VERSION})")
endif()
message(STATUS "Found VOTCA_XTP: ${VOTCA_XTP_DIR}")

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/xtp-tutorials
	PATTERN ".git*" EXCLUDE
	PATTERN "CMake*" EXCLUDE PATTERN "cmake_install.cmake" EXCLUDE
	PATTERN "install_manifest.txt" EXCLUDE PATTERN "*~" EXCLUDE)

configure_file(${PROJECT_SOURCE_DIR}/CMakeModules/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)

option(ENABLE_TESTING "Test if cmake worked" OFF)

if(ENABLE_TESTING)
  enable_testing()
  add_test(regression_xtp-tutorial_UninstallExists ${CMAKE_COMMAND} -DFileToCheck=${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
           -P ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/FileExists.cmake)

  set_tests_properties(regression_xtp-tutorial_UninstallExists PROPERTIES LABELS "xtp-tutorials")
endif(ENABLE_TESTING)
add_custom_target(uninstall-xtp-tutorials COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
if(NOT TARGET uninstall)
  add_custom_target(uninstall)
endif()
add_dependencies(uninstall uninstall-xtp-tutorials)

add_subdirectory(share/doc)
