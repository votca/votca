add_custom_target(gitversion-xtp COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)

file(GLOB_RECURSE VOTCA_SOURCES *.cc *.cpp)
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc )
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

if(NOT USE_GPU)
  list(FILTER VOTCA_SOURCES EXCLUDE REGEX ".*eigencuda.cc")
endif()

add_library(votca_xtp  ${VOTCA_SOURCES})
set_target_properties(votca_xtp PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(votca_xtp gitversion-xtp)
target_link_libraries(votca_xtp PUBLIC VOTCA::votca_csg VOTCA::votca_tools Boost::boost Eigen3::Eigen ${HDF5_LIBRARIES} PRIVATE LIBXC::LIBXC Boost::program_options Boost::serialization Boost::filesystem Boost::system Boost::timer )
if(MKL_FOUND)
  target_link_libraries(votca_xtp PRIVATE MKL::MKL)
endif()

if(USE_GPU)
  target_link_libraries(votca_xtp PRIVATE ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
endif()

#CMAKE_CURRENT_BINARY_DIR for gitversion.h
#PROJECT_BINARY_DIR/include for votca_config.h
target_include_directories(votca_xtp PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> ${HDF5_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS votca_xtp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
