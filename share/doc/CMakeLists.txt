if(NOT VOTCA_SPHINX_DIR)
  set(VOTCA_SPHINX_DIR "${CMAKE_CURRENT_BINARY_DIR}")
endif()

set(NOTEBOOKS "GROMACS/Methane/QMMM_GROMACS.ipynb")
list(APPEND NOTEBOOKS "LAMMPS/Thiophene/QMMM_LAMMPS.ipynb")
list(APPEND NOTEBOOKS "tools/dftgwbse_CH4/DFTGWBSE_ENERGY.ipynb")
list(APPEND NOTEBOOKS "tools/dftgwbse_CO_geoopt/DFTGWBSE_OPTIMIZATION.ipynb")
# list(APPEND NOTEBOOKS "GROMACS/KMC_Methane/GROMACS_KMC.ipynb")
list(APPEND NOTEBOOKS "LAMMPS/KMC_Thiophene/LAMMPS_KMC.ipynb")

find_package(JUPYTER_NBCONVERT 5.4.1)
find_package(SPHINX 2.1)
set(NEW_NOTEBOOKS)

if(JUPYTER_NBCONVERT_FOUND AND SPHINX_FOUND)
  enable_language(CXX) # to detect xtp
  find_package(VOTCA_XTP REQUIRED NO_MODULE)
  set_package_properties(VOTCA_XTP PROPERTIES TYPE REQUIRED PURPOSE "Votca xtp library")
  if(DEFINED VOTCA_XTP_VERSION AND NOT VOTCA_XTP_VERSION STREQUAL PROJECT_VERSION)
    message(FATAL_ERROR "Incompatible VOTCA XTP version found (needed ${PROJECT_VERSION}, found ${VOTCA_XTP_VERSION})")
  endif()
  message(STATUS "Found VOTCA_XTP: ${VOTCA_XTP_DIR}")

  get_target_property(XTP_PATH VOTCA::votca_xtp XTP_PATH)
  get_target_property(VOTCASHARE VOTCA::votca_xtp VOTCA_SHARE)
  foreach(_notebook ${NOTEBOOKS})
    get_filename_component(_FILE "${_notebook}" NAME_WE)
    get_filename_component(_DIR "${_notebook}" DIRECTORY)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_DIR})
    file(MAKE_DIRECTORY ${VOTCA_SPHINX_DIR}/xtp-tutorials)
    file(GLOB FILES ${PROJECT_SOURCE_DIR}/${_DIR}/*)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_notebook}
      COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/${_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${_DIR}
      DEPENDS ${FILES})
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_DIR}
      COMMAND ${CMAKE_COMMAND} -E env PATH=${XTP_PATH}:$ENV{PATH} VOTCASHARE=${VOTCASHARE} ${JUPYTER_EXECUTABLE} nbconvert --to notebook --ExecutePreprocessor.timeout=5000 --execute ${_FILE} --output ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_notebook}
    )
    list(APPEND NEW_NOTEBOOKS ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE})
  endforeach()
else()
  foreach(_DEP ${NOTEBOOKS})
    get_filename_component(_FILE "${_DEP}" NAME)
    add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_SOURCE_DIR}/${_DEP} ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE} DEPENDS ${PROJECT_SOURCE_DIR}/${_DEP}
)
    list(APPEND NEW_NOTEBOOKS ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE})
    endforeach()
endif()

add_custom_target(doc-xtp-tutorials DEPENDS ${NEW_NOTEBOOKS})

