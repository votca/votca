# Build step of BuildLibint.cmake, run by `cmake --build` (target votca_libint):
# compiles the extracted libint export and installs it into the per-machine
# cache. Holds the cache lock, so concurrent builds in other build
# directories wait and then find the finished install.
#   cmake -DSRC= -DPREFIX= -DSTAMP= -DLOCK= -DWORK= -DINIT= -DGENERATOR= -DJOBS= -DSHA= -P BuildLibintStep.cmake
file(LOCK "${LOCK}" GUARD PROCESS TIMEOUT 14400 RESULT_VARIABLE _lockres)
if(NOT _lockres EQUAL 0)
  message(FATAL_ERROR "Could not lock ${LOCK}: ${_lockres}")
endif()
if(EXISTS "${STAMP}")
  message(STATUS "libint2: built meanwhile by another build; using ${PREFIX}")
  return()
endif()

if(NOT EXISTS "${SRC}/CMakeLists.txt")
  message(FATAL_ERROR "libint2: extracted export ${SRC} is missing (cache removed?); rerun cmake")
endif()
message(STATUS "==== libint2: building into ${PREFIX} (once per machine; the VOTCA build continues afterwards) ====")
file(REMOVE_RECURSE "${WORK}")
file(GLOB _old LIST_DIRECTORIES true "${PREFIX}/*")
if(_old)
  file(REMOVE_RECURSE ${_old})
endif()

macro(_step)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE _rc)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "libint2: build step failed (${_rc}): ${ARGN}")
  endif()
endmacro()

_step(${CMAKE_COMMAND} -S "${SRC}" -B "${WORK}" -G "${GENERATOR}" -C "${INIT}"
  -DCMAKE_BUILD_TYPE=Release
  "-DCMAKE_INSTALL_PREFIX=${PREFIX}"
  -DCMAKE_INSTALL_LIBDIR=lib
  -DBUILD_SHARED_LIBS=OFF
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  -DLIBINT2_REQUIRE_CXX_API=ON
  -DLIBINT2_REQUIRE_CXX_API_COMPILED=OFF
  -DLIBINT2_ENABLE_FORTRAN=OFF
  -DLIBINT2_ENABLE_PYTHON=OFF
  -DBUILD_TESTING=OFF)
_step(${CMAKE_COMMAND} --build "${WORK}" --parallel ${JOBS} --config Release)
_step(${CMAKE_COMMAND} --install "${WORK}" --config Release)
file(WRITE "${STAMP}" "sha256=${SHA}\nsource=${SRC}\n")
file(REMOVE_RECURSE "${WORK}")
message(STATUS "==== libint2: installed into ${PREFIX}; continuing with VOTCA ====")
