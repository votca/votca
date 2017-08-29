cmake_minimum_required(VERSION 2.6)

project(csg-tutorials)

set(PROJECT_VERSION "1.5-dev")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

option(ENABLE_TESTING "Build and copy testing stuff" OFF)
if(VOTCA_PATH AND VOTCA_SHARE AND VOTCA_CSG_INVERSE AND ENABLE_TESTING)
  enable_testing()
  file(GLOB_RECURSE FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} settings.xml)
  foreach(DIR ${FILES})
    if(DIR MATCHES "(lammps|espresso|dlpoly|hoomd|atomistic|single_chain)")
      continue()
    endif()
    string(REGEX MATCH ^.*/ DIR ${DIR})
    list(APPEND DIRS ${DIR})
  endforeach(DIR)
  list(REMOVE_DUPLICATES DIRS)
  foreach(DIR ${DIRS})
    if (NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DIR})
      file(GLOB_RECURSE FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${DIR}/*)
      foreach(FILE ${FILES})
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/${FILE} ${CMAKE_CURRENT_BINARY_DIR}/${FILE})
      endforeach(FILE ${FILES})
    endif (NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    string(REGEX REPLACE "/$" "" DIR "${DIR}")
    string(REPLACE / _ TESTNAME "${DIR}")
    add_test(NAME ${TESTNAME} COMMAND ${VOTCA_CSG_INVERSE} --options settings.xml --do-iterations 1 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DIR})
    # run tests for csg (to test csg) and tools (coverage) as well
    set_tests_properties(${TESTNAME} PROPERTIES LABELS "csg-tutorials;csg;tools;votca")
    set(CSG_MDRUN_STEPS 500)
    if(DIR MATCHES "/(imc|realtime)$")
      set(CSG_MDRUN_STEPS 5000)
    elseif(DIR MATCHES "/cma/")
      set(CSG_MDRUN_STEPS 5000)
    elseif(DIR MATCHES "/re$")
      set(CSG_MDRUN_STEPS 2500)
    endif()
    set_tests_properties(${TESTNAME} PROPERTIES ENVIRONMENT "CSG_RUNTEST=yes;VOTCASHARE=${VOTCA_SHARE};CSG_MDRUN_STEPS=${CSG_MDRUN_STEPS};PATH=${VOTCA_PATH}:$ENV{PATH}")
  endforeach(DIR ${DIRS})
  message(STATUS "Use VOTCASHARE=${VOTCA_SHARE} PATH=${VOTCA_PATH} for out-of-cmake testing")
elseif(ENABLE_TESTING)
	message(FATAL_ERROR "We need VOTCA_PATH, VOTCA_SHARE and VOTCA_CSG_INVERSE set for testing!")
endif()

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ DESTINATION share/votca/csg-tutorials
	PATTERN ".git*" EXCLUDE PATTERN "step_*" EXCLUDE
	PATTERN "CMake*" EXCLUDE PATTERN "cmake_install.cmake" EXCLUDE
	PATTERN "install_manifest.txt" EXCLUDE PATTERN "*~" EXCLUDE
	PATTERN "inverse.log" EXCLUDE)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall-csg-tutorials COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
if(NOT TARGET uninstall)
  add_custom_target(uninstall)
endif()
add_dependencies(uninstall uninstall-csg-tutorials)
