install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/xtp-tutorials USE_SOURCE_PERMISSIONS
	PATTERN ".git*" EXCLUDE
	PATTERN "CMake*" EXCLUDE PATTERN "cmake_install.cmake" EXCLUDE
	PATTERN "install_manifest.txt" EXCLUDE PATTERN "*~" EXCLUDE)

if(NOT SPHINX_FOUND)
  return()
endif()
find_package(JUPYTER_NBCONVERT 5.4.1)
if(NOT JUPYTER_NBCONVERT_FOUND)
  return()
endif()

set(NOTEBOOKS "GROMACS/Methane/QMMM_GROMACS.ipynb")
list(APPEND NOTEBOOKS "LAMMPS/Thiophene/QMMM_LAMMPS.ipynb")
list(APPEND NOTEBOOKS "tools/dftgwbse_CH4/DFTGWBSE_ENERGY.ipynb")
list(APPEND NOTEBOOKS "LAMMPS/KMC_Thiophene/LAMMPS_KMC.ipynb")
list(APPEND NOTEBOOKS "pyxtp/geometry_optimization.ipynb")
list(APPEND NOTEBOOKS "pyxtp/single_point.ipynb")

set(NEW_NOTEBOOKS)

get_target_property(XTP_PATH VOTCA::votca_xtp XTP_PATH)
get_target_property(VOTCASHARE VOTCA::votca_xtp VOTCA_SHARE)
get_target_property(XTP_PYPATH VOTCA::xtp_binds XTP_PYPATH)
foreach(_notebook ${NOTEBOOKS})
  get_filename_component(_FILE "${_notebook}" NAME)
  get_filename_component(_DIR "${_notebook}" DIRECTORY)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_DIR})
  file(GLOB FILES ${CMAKE_CURRENT_SOURCE_DIR}/${_DIR}/*)
  add_custom_command(OUTPUT ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/${_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${_DIR}
    COMMAND ${CMAKE_COMMAND} -E env PATH=${XTP_PATH}:$ENV{PATH} VOTCASHARE=${VOTCASHARE} PYTHONPATH=${XTP_PYPATH}:$ENV{PYTHONPATH}  ${JUPYTER_EXECUTABLE} nbconvert --Application.log_level=0 --to notebook --ExecutePreprocessor.timeout=5000 --execute ${_FILE} --output ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_DIR}
    DEPENDS ${FILES}
  )
  list(APPEND NEW_NOTEBOOKS ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE})
endforeach()

add_custom_target(xtp-tutorials-rst DEPENDS ${NEW_NOTEBOOKS})
add_dependencies(sphinx-deps xtp-tutorials-rst)
