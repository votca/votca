cmake_minimum_required(VERSION 3.1)

project(votca-csg-manual NONE)

set(PROJECT_VERSION "1.5-dev")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

########################################################################
# User input options                                                   #
########################################################################
if (NOT DEFINED DATA)
  set(DATA "share/votca")
endif(NOT DEFINED DATA)

########################################################################
#Find external packages
########################################################################
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

find_package(TXT2TAGS REQUIRED)

find_program(CSG_CALL csg_call)
if(NOT CSG_CALL)
  message(FATAL_ERROR "csg_call not found")
endif()
find_program(CSG_PROPERTY csg_property)
if(NOT CSG_PROPERTY)
  message(FATAL_ERROR "csg_property not found")
endif()
execute_process(COMMAND ${CSG_CALL} --show-share OUTPUT_VARIABLE CSG_SHARE OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "CSG_SHARE is ${CSG_SHARE}")

set(TEX2EPS ${CMAKE_CURRENT_SOURCE_DIR}/scripts/tex2eps)
set(BBOX_ADD ${CMAKE_CURRENT_SOURCE_DIR}/scripts/bbox_add.pl)

include(${CMAKE_MODULE_PATH}/UseLATEX.cmake)

######################################
# Include the following subdirectory # 
######################################
#cannot use add_subdir as file depends only work with in the same directory
foreach(CURRENT_DIR functionality/fig usage/fig reference reference/xml)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR})
  include(${CURRENT_DIR}/CMakeLists.txt)
endforeach()

add_custom_target(gitversion COMMAND ${CMAKE_COMMAND}
  -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DPROJECT_VERSION="${PROJECT_VERSION}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitid.tex)
add_custom_command(OUTPUT gitid.tex DEPENDS gitversion)

configure_file(${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
add_custom_target(uninstall-csg-manual COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
if(NOT TARGET uninstall)
  add_custom_target(uninstall)
endif()
add_dependencies(uninstall uninstall-csg-manual)

set(LATEX_DEFAULT_BUILD ps)
add_latex_document(manual.tex
  BIBFILES votca.bib
  INPUTS inputfiles.tex introduction.tex newcommands.tex reference.tex theory.tex titlepage.tex usage.tex
    usage/boltzmann.tex usage/cgruns.tex usage/cibi.tex usage/dlpoly.tex usage/forcematching.tex usage/iterative_methods.tex
    functionality/propane.xml functionality/settings.xml
  IMAGES fig/logo.eps fig/propane.eps usage/fig/propane.eps fig/excl.eps fig/hexane2.eps fig/propane_hist2d.eps
  TARGET_NAME manual DEPENDS gitid.tex functionality/fig/flowchart.eps functionality/fig/iterativemethods.eps
    usage/fig/flow_boltzmann.eps usage/fig/flow_fmatch.eps usage/fig/flow_ibi.eps reference/csgid.tex
    reference/xml/mapping.xml.tex reference/xml/csg_defaults.xml.tex reference/xml/topol.xml.tex
)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  include(FeatureSummary)
  feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
endif (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
