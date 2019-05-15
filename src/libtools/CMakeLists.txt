find_package(EXPAT REQUIRED)

file(GLOB VOTCA_SOURCES *.cc)
file(GLOB NOT_VOTCA_SOURCES version_nb.cc)
file(GLOB_RECURSE VOTCA_SQL_SOURCES database.cc statement.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES} ${VOTCA_SQL_SOURCES})

add_library(votca_tools ${VOTCA_SOURCES} ${VOTCA_LINALG_SOURCES})

option(WITH_SQLITE3 "Use SQLite3, at this point only needed for charge transport parts!" ON)
if (WITH_SQLITE3)
  find_package(SQLITE3 REQUIRED)
  target_sources(votca_tools PRIVATE ${VOTCA_SQL_SOURCES})
  target_link_libraries(votca_tools PRIVATE SQLite::SQLite3)
endif(WITH_SQLITE3)


#for gitversion.h
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(gitversion COMMAND ${CMAKE_COMMAND}
  -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -DMERCURIAL_EXECUTABLE="${MERCURIAL_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)

add_dependencies(votca_tools gitversion)

set_target_properties(votca_tools PROPERTIES SOVERSION ${SOVERSION})
target_link_libraries(votca_tools PUBLIC Boost::boost Boost::program_options Eigen3::Eigen PRIVATE
	Boost::filesystem Boost::system
	FFTW3::FFTW3 EXPAT::EXPAT Threads::Threads ${MATH_LIBRARIES})
if(MKL_FOUND)
  target_link_libraries(votca_tools PRIVATE MKL::MKL)
endif()
install(TARGETS votca_tools LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
