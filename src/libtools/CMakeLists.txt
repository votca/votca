find_package(EXPAT)
if(NOT EXPAT_FOUND)
  message(FATAL_ERROR "EXPAT not found,  make sure you have also installed the expat and it's dev package")
endif(NOT EXPAT_FOUND)
include_directories(${EXPAT_INCLUDE_DIR})
set (EXPAT_CFLAGS_PKG "-I${EXPAT_INCLUDE_DIRS}")
string(REGEX REPLACE "^(.*)/lib([^/]*)\\.[^.]*$" "-L\\1 -l\\2" EXPAT_LIBS_PKG ${EXPAT_LIBRARY})

option(WITH_FFTW "Use FFTW3, disabling leads to reduced functionality!" ON)
if (WITH_FFTW)
  find_package(FFTW3)
  if(NOT FFTW3_FOUND)
    message(FATAL_ERROR "FFTW3 not found,  make sure you have also installed the fftw3 and it's dev package")
  endif(NOT FFTW3_FOUND)
  include_directories(${FFTW3_INCLUDE_DIRS})
else(WITH_FFTW)
  #used in votca_config.h
  set(NOFFTW TRUE)
endif(WITH_FFTW)

option(WITH_GSL "Use GSL, disabling leads to reduced functionality!" ON)
if (WITH_GSL)
  find_package(GSL)
  if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL not found,  make sure you have also installed the gsl and it's dev package")
  endif(NOT GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
else(WITH_GSL)
  #used in votca_config.h
  set(NOGSL TRUE)
endif(WITH_GSL)

configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(hgversion COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${CMAKE_MODULE_PATH}/hgversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES hgversion.h)

file(GLOB_RECURSE VOTCA_SOURCES *.cc)
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc test.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})
add_library(votca_tools ${VOTCA_SOURCES})
add_dependencies(votca_tools hgversion)
target_link_libraries(votca_tools ${BOOST_LIBRARIES} ${GSL_LIBRARIES}
  ${FFTW3_LIBRARIES} ${EXPAT_LIBRARIES} ${THREAD_LIBRARIES} ${MATH_LIBRARIES})
install(TARGETS votca_tools LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_tools.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_tools.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_tools.pc DESTINATION ${LIB}/pkgconfig)

add_executable(linktest linktest.cc)
target_link_libraries(linktest votca_tools)
