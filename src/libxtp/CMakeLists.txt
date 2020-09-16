add_custom_target(gitversion-xtp COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${PROJECT_SOURCE_DIR}/CMakeModules/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)

file(GLOB_RECURSE VOTCA_SOURCES *.cc *.cpp)
file(GLOB_RECURSE NOT_VOTCA_SOURCES cudapipeline.cc cudamatrix.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

get_filename_component(_PARENT_PATH "${CMAKE_CURRENT_BINARY_DIR}" DIRECTORY)

add_library(votca_xtp  ${VOTCA_SOURCES})
set_target_properties(votca_xtp PROPERTIES SOVERSION ${SOVERSION})
set_target_properties(votca_xtp PROPERTIES VOTCA_SHARE ${PROJECT_SOURCE_DIR}/share)
set_target_properties(votca_xtp PROPERTIES XTP_PATH ${_PARENT_PATH}/tools)
add_dependencies(votca_xtp gitversion-xtp)
target_link_libraries(votca_xtp PUBLIC VOTCA::votca_csg VOTCA::votca_tools Boost::boost Eigen3::Eigen ${HDF5_LIBRARIES} PRIVATE PkgConfig::LIBXC PkgConfig::LIBINT Boost::program_options Boost::filesystem Boost::system Boost::timer )

if(OpenMP_FOUND)
  target_link_libraries(votca_xtp PUBLIC OpenMP::OpenMP_CXX)
endif()

if(USE_CUDA)
  target_link_libraries(votca_xtp PUBLIC ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
  target_include_directories(votca_xtp PUBLIC ${CUDA_INCLUDE_DIRS})
  target_sources(votca_xtp PRIVATE cudapipeline.cc cudamatrix.cc)
endif()

#CMAKE_CURRENT_BINARY_DIR for gitversion.h
#PROJECT_BINARY_DIR/include/votca/xtp for votca_xtp_config.h
target_include_directories(votca_xtp PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> ${HDF5_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/votca/xtp>
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
