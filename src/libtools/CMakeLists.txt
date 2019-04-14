find_package(EXPAT)
if(NOT EXPAT_FOUND)
  message(FATAL_ERROR "EXPAT not found,  make sure you have also installed the expat and it's dev package")
endif(NOT EXPAT_FOUND)
include_directories(${EXPAT_INCLUDE_DIR})
set (EXPAT_CFLAGS_PKG "-I${EXPAT_INCLUDE_DIRS}")
string(REGEX REPLACE "^(.*)/lib([^/]*)\\.[^.]*$" "-L\\1 -l\\2" EXPAT_LIBS_PKG ${EXPAT_LIBRARY})

file(GLOB VOTCA_SOURCES *.cc)
file(GLOB NOT_VOTCA_SOURCES version_nb.cc)
file(GLOB_RECURSE VOTCA_SQL_SOURCES database.cc statement.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES} ${VOTCA_SQL_SOURCES})

#for gitversion.h
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(gitversion COMMAND ${CMAKE_COMMAND}
  -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -DMERCURIAL_EXECUTABLE="${MERCURIAL_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)

add_library(votca_tools ${VOTCA_SOURCES} ${VOTCA_SQL_SOURCES} ${VOTCA_LINALG_SOURCES})
add_dependencies(votca_tools gitversion)

set_target_properties(votca_tools PROPERTIES SOVERSION ${SOVERSION})
target_link_libraries(votca_tools ${Boost_LIBRARIES} ${MKL_LIBRARIES} ${SQLITE3_LIBRARIES}
  ${FFTW3_LIBRARIES} ${EXPAT_LIBRARIES} ${THREAD_LIBRARIES} ${MATH_LIBRARIES})
install(TARGETS votca_tools LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

configure_file(libvotca_tools.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_tools.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_tools.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
