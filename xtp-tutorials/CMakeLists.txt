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

set(NEW_NOTEBOOKS)

# Run the notebooks ONE AT A TIME, not as many at once as the build's -j
# allows.
#
# Every notebook below executes real calculations -- full DFT, GW-BSE,
# QMMM or KMC -- and each of those helps itself to threads twice over:
# xtp sizes its OpenMP pool from the machine's core count (CI logs show
# "Using 4 threads" and "Davidson Solver using 4 threads"), and several
# cells add process-level parallelism on top with xtp_parallel -x 2. The
# CI runner then sets -j to ncpu+1 (.github/actions/setup/setup.sh), so
# Ninja is free to start five notebooks at once on a four-core runner:
# roughly forty threads fighting over four cores, plus the memory of five
# concurrent GW-BSE runs.
#
# That is the shape of a flake that "usually passes on a re-run": nothing
# is wrong with the notebooks, they are simply being starved, and the
# first one to exceed --ExecutePreprocessor.timeout=5000 takes the whole
# sphinx target down with it. Interleaved output from five notebooks also
# makes the failure nearly impossible to read afterwards -- which is its
# own cost, separate from the failure itself.
#
# A job pool of one serializes just these commands while leaving the rest
# of the build's parallelism untouched. Each notebook then gets the whole
# machine, so serializing them costs far less wall-clock than the ncpu+1
# arithmetic suggests. JOB_POOL is honoured by Ninja (the generator CI
# uses) and silently ignored by the others, so this is safe either way.
set_property(GLOBAL APPEND PROPERTY JOB_POOLS votca_notebooks=1)

get_target_property(XTP_PATH VOTCA::votca_xtp XTP_PATH)
get_target_property(VOTCASHARE VOTCA::votca_xtp VOTCA_SHARE)

if(TARGET VOTCA::xtp_binds)
  get_target_property(XTP_PYPATH VOTCA::xtp_binds XTP_PYPATH)
  list(APPEND NOTEBOOKS "pyxtp/geometry_optimization.ipynb")
  list(APPEND NOTEBOOKS "pyxtp/single_point.ipynb")
endif()

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
    JOB_POOL votca_notebooks
  )
  list(APPEND NEW_NOTEBOOKS ${VOTCA_SPHINX_DIR}/xtp-tutorials/${_FILE})
endforeach()

add_custom_target(xtp-tutorials-rst DEPENDS ${NEW_NOTEBOOKS})
add_dependencies(sphinx-deps xtp-tutorials-rst)
