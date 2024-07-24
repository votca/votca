file(GLOB_RECURSE VOTCA_SOURCES *.cc *.cpp)
file(GLOB_RECURSE NOT_VOTCA_SOURCES cudapipeline.cc cudamatrix.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

get_filename_component(_PARENT_PATH "${CMAKE_CURRENT_BINARY_DIR}" DIRECTORY)

add_library(votca_xtp  ${VOTCA_SOURCES})
target_compile_features(votca_xtp PUBLIC cxx_std_17)

set_target_properties(votca_xtp PROPERTIES SOVERSION ${SOVERSION})
set_target_properties(votca_xtp PROPERTIES VOTCA_SHARE ${PROJECT_SOURCE_DIR}/xtp/share)
set_target_properties(votca_xtp PROPERTIES XTP_PATH ${_PARENT_PATH}/tools)
target_link_libraries(votca_xtp PUBLIC VOTCA::votca_csg VOTCA::votca_tools Boost::boost Eigen3::Eigen ${HDF5_LIBRARIES} PRIVATE Libxc::xc Libint2::int2 ECPINT::ecpint Boost::program_options Boost::filesystem Boost::system Boost::timer Boost::serialization)

if(USE_CUDA)
  target_link_libraries(votca_xtp PUBLIC ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
  target_include_directories(votca_xtp PUBLIC ${CUDA_INCLUDE_DIRS})
  target_sources(votca_xtp PRIVATE cudapipeline.cc cudamatrix.cc)
endif()

# PROJECT_BINARY_DIR/xtp/include/votca/xtp for votca_xtp_config.h
target_include_directories(votca_xtp PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/xtp/include> ${HDF5_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/xtp/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/xtp/include/votca/xtp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS votca_xtp EXPORT VOTCA_XTP_Targets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT VOTCA_XTP_Targets FILE VOTCA_XTP_Targets.cmake NAMESPACE VOTCA:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VOTCA_XTP)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("VOTCA_XTPConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
configure_file(VOTCA_XTPConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_XTPConfig.cmake" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_XTPConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_XTPConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VOTCA_XTP)
add_library(VOTCA::votca_xtp ALIAS votca_xtp)
