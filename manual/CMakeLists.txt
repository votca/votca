cmake_minimum_required(VERSION 3.1)

project(votca-xtp-manual NONE)

set(PROJECT_VERSION "1.5-dev")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

########################################################################
# User input options                                                   #
########################################################################
set(CMAKE_SIZEOF_VOID_P 8)
include(GNUInstallDirs)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
	message(FATAL_ERROR "In-source build not supported with UseLATEX.cmake\nDo 'rm CMakeCache.txt; mkdir build; cd build; cmake ..' instead")
endif()

########################################################################
#Find external packages
########################################################################
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)

find_package(TXT2TAGS REQUIRED)

foreach(prog csg_property dvipdf inkscape votca_property)
  string(TOUPPER "${prog}" PROG)
  if(TARGET ${prog})
    set(${PROG} $<TARGET_FILE:${prog}>)
  endif()
  find_program(${PROG} ${prog})
  find_package_handle_standard_args(${PROG} REQUIRED_VARS ${PROG})
  if(NOT ${PROG}_FOUND)
    message(FATAL_ERROR "Could not find ${prog}")
  endif()
endforeach()

set(FIGTEX2EPS ${CMAKE_CURRENT_SOURCE_DIR}/scripts/figtex2eps)

include(${CMAKE_MODULE_PATH}/UseLATEX.cmake)

######################################
# Include the following subdirectory # 
######################################
#cannot use add_subdir as file depends only work with in the same directory
foreach(CURRENT_DIR fig/chemical_structure fig/idft_flow fig/stochastic
  fig/reorganization_energy fig/workflow reference/xml reference/programs
  reference/calculators)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR})
  include(${CURRENT_DIR}/CMakeLists.txt)
endforeach()

add_custom_target(xtp-gitid COMMAND ${CMAKE_COMMAND}
  -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DPROJECT_VERSION="${PROJECT_VERSION}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitid.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitid.tex)
add_custom_command(OUTPUT gitid.tex DEPENDS xtp-gitid)

add_latex_document(xtp-manual.tex
  BIBFILES xtp-manual.bib FORCE_DVI
  INPUTS newcommands.tex titlepage.tex theory/introduction.tex theory/workflow.tex theory/morphology.tex
    theory/fragments.tex theory/neighborlist.tex theory/reorganization_energy.tex theory/energetic_disorder.tex
    theory/gdma.tex theory/thole.tex theory/electronic_couplings.tex theory/dipro.tex theory/dft.tex
    theory/moo.tex theory/rate.tex theory/kmc.tex theory/nondispersive.tex theory/stochastic.tex
    theory/analysis.tex theory/correlations.tex input/topology_atomistic.tex input/topology_coarsegrained.tex
    input/orbitals_zindo.tex input/header_edft.tex input/header_idft.tex input/output_dft.tex
    input/statefile.tex reference/reference.tex reference/calculators.xml input/dcv2t/dcv2t.pdb
    input/dcv2t/map.xml input/get_orbitals.com input/qmpackage_headers/gaussian_edft.xml
    input/qmpackage_headers/turbomole_edft.xml input/qmpackage_headers/nwchem_edft.xml input/TI.xml
    input/qmpackage_headers/gaussian_idft.xml input/qmpackage_headers/turbomole_idft.xml
    input/qmpackage_headers/nwchem_idft.xml fig/stochastic/stochastic.eps
    fig/conjugated_segment/fragment_segment.eps RandomFacts.tex
  IMAGES fig/logo.eps fig/logo_trans.eps
  TARGET_NAME xtp-manual DEPENDS gitid.tex fig/chemical_structure/dcv2t.eps fig/idft_flow/scheme.eps
    fig/stochastic/overview_stochastic.eps fig/idft_flow/schemes_all.eps
    fig/reorganization_energy/marcus_parabolas.eps fig/reorganization_energy/monomer_parabolas.eps
    fig/workflow/workflow.eps reference/xml/map.xml.tex reference/xml/options.xml.tex
    reference/programs/all.tex reference/calculators/calculators.tex
)

add_custom_command(OUTPUT xtp-manual.pdf
  COMMAND ${DVIPDF} xtp-manual.dvi xtp-manual.pdf
  DEPENDS xtp-manual_dvi
)
add_custom_target(xtp-manual_pdf ALL DEPENDS xtp-manual.pdf)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/xtp-manual.pdf DESTINATION ${CMAKE_INSTALL_DOCDIR})
