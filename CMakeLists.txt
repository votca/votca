cmake_minimum_required(VERSION 2.6)

project(csg-tutorials)

set(PROJECT_VERSION "1.3-dev")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ DESTINATION share/votca/csg-tutorials
	PATTERN ".hg*" EXCLUDE PATTERN "step_*" EXCLUDE
	PATTERN "CMake*" EXCLUDE PATTERN "cmake_install.cmake" EXCLUDE
	PATTERN "install_manifest.txt" EXCLUDE PATTERN "*~" EXCLUDE
	PATTERN "inverse.log" EXCLUDE)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
