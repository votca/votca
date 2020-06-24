if(NOT VOTCA_SPHINX_DIR)
  set(VOTCA_SPHINX_DIR "${CMAKE_CURRENT_BINARY_DIR}")
endif()

file(GLOB NOTEBOOKS "${PROJECT_SOURCE_DIR}/GROMACS/Methane/QMMM_GROMACS.ipynb")
list(APPEND NOTEBOOKS "${PROJECT_SOURCE_DIR}/LAMMPS/Thiophene/QMMM_LAMMPS.ipynb")
list(APPEND NOTEBOOKS "${PROJECT_SOURCE_DIR}/tools/dftgwbse_CH4/DFTGWBSE_ENERGY.ipynb")
list(APPEND NOTEBOOKS "${PROJECT_SOURCE_DIR}/tools/dftgwbse_CO_geoopt/DFTGWBSE_OPTIMIZATION.ipynb")
# list(APPEND NOTEBOOKS "${PROJECT_SOURCE_DIR}/GROMACS/KMC_Methane/GROMACS_KMC.ipynb")
list(APPEND NOTEBOOKS "${PROJECT_SOURCE_DIR}/LAMMPS/KMC_Thiophene/LAMMPS_KMC.ipynb")

find_program(JUPYTER_EXECUTABLE NAMES jupyter "Interactive computing environment (https://jupyter.org/)")
find_program(PIP_EXECUTABLE NAMES pip "Python package installer (https://pip.pypa.io/en/stable/)")

add_custom_target(jupyter_client COMMAND ${PIP_EXECUTABLE} install wheel ipykernel jupyter_client --user)

get_target_property(XTP_PATH votca_xtp XTP_PATH)
get_target_property(VOTCASHARE votca_xtp VOTCA_SHARE)

set(NEW_NOTEBOOKS)
foreach(_notebook ${NOTEBOOKS})
  get_filename_component(_FILE "${_notebook}" NAME_WE)
  get_filename_component(_DIR "${_notebook}" DIRECTORY)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/${_FILE}
    WORKING_DIRECTORY ${_DIR}
    COMMAND ${CMAKE_COMMAND} -E env PATH=${XTP_PATH}:$ENV{PATH} VOTCASHARE=${VOTCASHARE} ${JUPYTER_EXECUTABLE} nbconvert --to notebook --ExecutePreprocessor.timeout=1000 --execute ${_notebook} --output ${VOTCA_SPHINX_DIR}/${_FILE}
    )
    list(APPEND NEW_NOTEBOOKS ${VOTCA_SPHINX_DIR}/${_FILE})
  endforeach()

add_custom_target(doc-xtp-tutorials DEPENDS ${NEW_NOTEBOOKS})
add_dependencies(doc-xtp-tutorials jupyter_client)
