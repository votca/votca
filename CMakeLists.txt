cmake_minimum_required(VERSION 3.12)

project(votca-csg-manual LANGUAGES NONE)

set(PROJECT_VERSION "2021-dev")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

########################################################################
# User input options                                                   #
########################################################################
set(CMAKE_SIZEOF_VOID_P 8)
include(GNUInstallDirs)
include(FeatureSummary)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source build not supported with UseLATEX.cmake\nDo 'rm CMakeCache.txt; mkdir build; cd build; cmake ..' instead")
endif()

########################################################################
#Find external packages
########################################################################
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
  set_package_properties(Git PROPERTIES TYPE OPTIONAL PURPOSE "Generated version for development version")
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

find_package(TXT2TAGS REQUIRED)
set_package_properties(TXT2TAGS PROPERTIES TYPE OPTIONAL PURPOSE "Used to build csg-manual")
find_package(VOTCA_CSG NO_MODULE REQUIRED)
set_package_properties(VOTCA_CSG PROPERTIES TYPE REQUIRED PURPOSE "Votca csg library")
if(DEFINED VOTCA_CSG_VERSION AND NOT VOTCA_CSG_VERSION STREQUAL PROJECT_VERSION)
  message(FATAL_ERROR "Incompatible VOTCA CSG version found (needed ${PROJECT_VERSION}, found ${VOTCA_CSG_VERSION})")
endif()
message(STATUS "Found VOTCA_CSG: ${VOTCA_CSG_DIR}")

get_target_property(CSG_SHARE VOTCA::csg_call VOTCASHARE)
message(STATUS "CSG_SHARE is ${CSG_SHARE}")

get_target_property(VOTCA_CSG_DEFAULTS VOTCA::csg_inverse VOTCA_CSG_DEFAULTS)
message(STATUS "VOTCA_CSG_DEFAULTS is ${VOTCA_CSG_DEFAULTS}")

add_custom_target(gitid COMMAND ${CMAKE_COMMAND}
  -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DPROJECT_VERSION="${PROJECT_VERSION}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${PROJECT_SOURCE_DIR}/CMakeModules/gitid.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitid.tex)
add_custom_command(OUTPUT gitid.tex DEPENDS gitid)

configure_file(${PROJECT_SOURCE_DIR}/CMakeModules/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall-csg-manual COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
if(NOT TARGET uninstall)
  add_custom_target(uninstall)
endif()
add_dependencies(uninstall uninstall-csg-manual)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/index.rst DESTINATION ${CMAKE_INSTALL_DOCDIR})

option(ENABLE_TESTING "Test if manual exists" OFF)
add_feature_info(ENABLE_TESTING ENABLE_TESTING "Enable unit tests")
if(ENABLE_TESTING)
  enable_testing()
  add_test(ManualExists ${CMAKE_COMMAND} -DFileToCheck=${CMAKE_CURRENT_SOURCE_DIR}/index.rst
           -P ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/FileExists.cmake) 
  set_tests_properties(ManualExists PROPERTIES LABELS "csg-manual;votca")
endif(ENABLE_TESTING)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
endif (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)

if(NOT VOTCA_SPHINX_DIR)
  set(VOTCA_SPHINX_DIR "${CMAKE_CURRENT_BINARY_DIR}")
endif()

file(GLOB RSTS ${CMAKE_CURRENT_SOURCE_DIR}/*.rst)
file(GLOB PNGS ${CMAKE_CURRENT_SOURCE_DIR}/fig/*.png)
set(MANUAL_RST_FILES)
set(MANUAL_PNG_FILES)
foreach(_RST ${RSTS})
  get_filename_component(_FILE "${_RST}" NAME)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg-manual/${_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_RST} ${VOTCA_SPHINX_DIR}/csg-manual/${_FILE}
    DEPENDS ${_RST})
  list(APPEND MANUAL_RST_FILES ${VOTCA_SPHINX_DIR}/csg-manual/${_FILE})
endforeach()
foreach(_PNG ${PNGS})
    get_filename_component(_FILE "${_PNG}" NAME)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/csg-manual/fig/${_FILE}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_PNG} ${VOTCA_SPHINX_DIR}/csg-manual/fig/${_FILE}
      DEPENDS ${_PNG})
  list(APPEND MANUAL_PNG_FILES ${VOTCA_SPHINX_DIR}/csg-manual/fig/${_FILE})
endforeach()
add_custom_target(doc-csg-manual DEPENDS ${MANUAL_RST_FILES} ${MANUAL_PNG_FILES})
