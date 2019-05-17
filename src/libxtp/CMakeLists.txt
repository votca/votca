add_custom_target(gitversion-xtp COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

file(GLOB_RECURSE VOTCA_SOURCES *.cc *.cpp)
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc )
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

add_library(votca_xtp  ${VOTCA_SOURCES})
set_target_properties(votca_xtp PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(votca_xtp gitversion-xtp)
target_link_libraries(votca_xtp PUBLIC VOTCA::votca_ctp Boost::boost PRIVATE ${MKL_LIBRARIES} LIBXC::LIBXC VOTCA::votca_csg VOTCA::votca_tools Boost::program_options Boost::filesystem Boost::system Boost::timer ${HDF5_LIBRARIES} )
install(TARGETS votca_xtp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
