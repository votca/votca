cmake_minimum_required(VERSION 3.10)

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

set(TEX2EPS ${CMAKE_CURRENT_SOURCE_DIR}/scripts/tex2eps)
set(BBOX_ADD ${CMAKE_CURRENT_SOURCE_DIR}/scripts/bbox_add.pl)

include(${PROJECT_SOURCE_DIR}/CMakeModules/UseLATEX.cmake)

if(NOT LATEX_WORKS)
  #test if latex actually works, this also generates latex.fmt on fresh docker containers
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/latex-tmp)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/latex-tmp/test.tex "\\documentclass{article}\n\\begin{document}\ntest\n\\end{document}\n")
  execute_process(COMMAND ${LATEX_COMPILER} -halt-on-error test.tex
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/latex-tmp
    OUTPUT_VARIABLE LATEX_OUTPUT ERROR_VARIABLE LATEX_OUTPUT
    RESULT_VARIABLE LATEX_EXIT)
  if(LATEX_EXIT EQUAL 0)
    message(STATUS "LaTeX worked")
    set(LATEX_WORKS TRUE CACHE INTERNAL "LaTeX works")
  else()
    message(FATAL_ERROR "LaTeX doesn't work, see ${CMAKE_CURRENT_BINARY_DIR}/latex-tmp for details")
  endif()
endif()

######################################
# Include the following subdirectory # 
######################################
#cannot use add_subdir as file depends only work with in the same directory
foreach(CURRENT_DIR functionality/fig usage/fig reference reference/xml
  reference/programs reference/scripts)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR})
  include(${CURRENT_DIR}/CMakeLists.txt)
endforeach()

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

add_latex_document(manual.tex
  BIBFILES votca.bib FORCE_DVI
  INPUTS inputfiles.tex introduction.tex newcommands.tex reference.tex theory.tex titlepage.tex usage.tex
    usage/boltzmann.tex usage/cgruns.tex usage/cibi.tex usage/dlpoly.tex usage/forcematching.tex usage/iterative_methods.tex
    functionality/propane.xml functionality/settings.xml advanced_topics.tex
  IMAGES fig/logo.eps fig/propane.eps usage/fig/propane.eps fig/excl.eps fig/hexane2.eps fig/propane_hist2d.eps
  TARGET_NAME csg-manual DEPENDS gitid.tex functionality/fig/flowchart.eps functionality/fig/iterativemethods.eps
    usage/fig/flow_boltzmann.eps usage/fig/flow_fmatch.eps usage/fig/flow_ibi.eps reference/csgid.tex
    reference/xml/mapping.xml.tex reference/xml/csg_defaults.xml.tex reference/xml/topol.xml.tex
    reference/programs/all.tex reference/scripts/csg_table.tex reference/scripts/all.tex
)

add_custom_command(OUTPUT manual.pdf
  COMMAND ${DVIPS_CONVERTER} manual.dvi
  COMMAND ${PS2PDF_CONVERTER} manual.ps manual.pdf
  DEPENDS csg-manual_dvi
)
add_custom_target(csg-manual_pdf ALL DEPENDS manual.pdf)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/manual.pdf DESTINATION ${CMAKE_INSTALL_DOCDIR})

option(ENABLE_TESTING "Test if manual exists" OFF)
add_feature_info(ENABLE_TESTING ENABLE_TESTING "Enable unit tests")
if(ENABLE_TESTING)
  enable_testing()
  add_test(ManualExists ${CMAKE_COMMAND} -DFileToCheck=${CMAKE_CURRENT_BINARY_DIR}/manual.pdf
           -P ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/FileExists.cmake) 
  set_tests_properties(ManualExists PROPERTIES LABELS "csg-manual;votca")
endif(ENABLE_TESTING)

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
endif (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
