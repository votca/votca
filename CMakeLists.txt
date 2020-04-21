cmake_minimum_required(VERSION 3.1)

project(xtp-tutorials)

set(PROJECT_VERSION "1.7-dev")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)
include(GNUInstallDirs)

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
