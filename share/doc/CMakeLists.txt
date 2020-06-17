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
set(NEW_NOTEBOOKS)
foreach(_notebook ${NOTEBOOKS})
  get_filename_component(_FILE "${_notebook}" NAME_WLE)
  get_filename_component(_DIR "${_notebook}" DIRECTORY)
  add_custom_command(OUTPUT ${_notebook}.ipynb
    COMMAND ${JUPYTER_EXECUTABLE} nbconvert --to notebook --ExecutePreprocessor.timeout=1000 --inplace --execute ${_notebook}
    )
    list(APPEND NEW_NOTEBOOKS ${_notebook}.ipynb)
endforeach()
add_custom_target(executed_notebooks DEPENDS ${NEW_NOTEBOOKS})


set(IPYNB_FILES)
foreach(_DEP ${NEW_NOTEBOOKS})  
  get_filename_component(_FILE "${_DEP}" NAME)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/${_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy ${_DEP} ${VOTCA_SPHINX_DIR}/${_FILE}
    DEPENDS executed_notebooks ${NEW_NOTEBOOKS}
    )
  list(APPEND IPYNB_FILES ${VOTCA_SPHINX_DIR}/${_FILE})
endforeach()

add_custom_target(doc-xtp-tutorials DEPENDS ${IPYNB_FILES})
