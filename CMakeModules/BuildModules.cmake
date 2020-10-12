include(ExternalProject)
set(ENABLED_VOTCA_PACKAGES tools csg csg-tutorials)
if(BUILD_CSG_MANUAL)
  list(APPEND ENABLED_VOTCA_PACKAGES csg-manual)
endif()
if(BUILD_CSGAPPS)
  list(APPEND ENABLED_VOTCA_PACKAGES csgapps)
endif()
if(BUILD_XTP)
  list(APPEND ENABLED_VOTCA_PACKAGES xtp xtp-tutorials)
endif()

include(GNUInstallDirs)
set(VOTCA_RPATH_OPTS)
# for manuals/csg-tutorials we need injected rpath to call binaries
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(VOTCA_RPATH_OPTS "-DCMAKE_INSTALL_RPATH=${CMAKE_INSTALL_FULL_LIBDIR}")
endif("${isSystemDir}" STREQUAL "-1")

message(STATUS "Building ${ENABLED_VOTCA_PACKAGES} as modules")
foreach(VOTCA_MODULE ${ENABLED_VOTCA_PACKAGES})
  if(VOTCA_TARBALL_DIR AND VOTCA_TARBALL_TAG)
    set(VOTCA_SOURCE_DIR URL ${VOTCA_TARBALL_DIR}/votca-${VOTCA_MODULE}-${VOTCA_TARBALL_TAG}.tar.gz)
  else()
    set(VOTCA_SOURCE_DIR SOURCE_DIR ${PROJECT_SOURCE_DIR}/${VOTCA_MODULE})
  endif()

  set(VOTCA_DEPENDS)
  set(VOTCA_EXTRA_OPTS -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS})
  set(VOTCA_TEST_COMMAND make test CTEST_OUTPUT_ON_FAILURE=1)
  if(VOTCA_MODULE STREQUAL "csg")
    set(VOTCA_DEPENDS tools)
  elseif(VOTCA_MODULE STREQUAL "csgapps")
    set(VOTCA_DEPENDS csg)
  elseif(VOTCA_MODULE STREQUAL "csg-manual")
    set(VOTCA_DEPENDS csg)
  elseif(VOTCA_MODULE STREQUAL "csg-tutorials")
    set(VOTCA_DEPENDS csg)
    list(REMOVE_ITEM VOTCA_EXTRA_OPTS -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS})
    if(VOTCA_TEST_OPTS)
      message(STATUS "Found VOTCA_TEST_OPTS: ${VOTCA_TEST_OPTS}")
      list(APPEND VOTCA_TEST_COMMAND "ARGS=${VOTCA_TEST_OPTS}")
    endif()
    if(ENABLE_REGRESSION_TESTING)
      list(APPEND VOTCA_EXTRA_OPTS -DENABLE_REGRESSION_TESTING:BOOL=${ENABLE_REGRESSION_TESTING})
    endif()
  elseif(VOTCA_MODULE STREQUAL "xtp")
    set(VOTCA_DEPENDS csg)
  elseif(VOTCA_MODULE STREQUAL "xtp-tutorials")
    set(VOTCA_DEPENDS xtp)
  endif()
  if(NOT ENABLE_TESTING)
    set(VOTCA_TEST_COMMAND)
  endif()
  ExternalProject_Add(${VOTCA_MODULE}
   DEPENDS ${VOTCA_DEPENDS}
   TEST_COMMAND ${VOTCA_TEST_COMMAND}
   TEST_BEFORE_INSTALL ${VOTCA_TEST_ENABLED}
   ${VOTCA_SOURCE_DIR}
   CMAKE_ARGS 
     -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
     -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
     -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
     -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
     -DCMAKE_PREFIX_PATH:PATH=${CMAKE_INSTALL_PREFIX}
     -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=${CMAKE_INSTALL_RPATH_USE_LINK_PATH}
     -DENABLE_TESTING:BOOL=${ENABLE_TESTING}
     ${VOTCA_RPATH_OPTS}
     ${VOTCA_EXTRA_OPTS}
  )
endforeach()

if(ENABLE_TESTING)
  add_test(IncludeExists ${CMAKE_COMMAND} -DFileToCheck=${CMAKE_INSTALL_FULL_INCLUDEDIR}/votca
           -P ${PROJECT_SOURCE_DIR}/CMakeModules/FileExists.cmake) 
  set_tests_properties(IncludeExists PROPERTIES LABELS "votca")
endif(ENABLE_TESTING)
install(FILES ${PROJECT_SOURCE_DIR}/README.rst DESTINATION ${CMAKE_INSTALL_DOCDIR})