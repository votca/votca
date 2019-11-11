file(GLOB GMX_SOURCES modules/io/gmx*.cc)

file(GLOB H5MD_SOURCES modules/io/h5md*.cc)

file(GLOB IO_SOURCES modules/io/*.cc)
list(REMOVE_ITEM IO_SOURCES ${GMX_SOURCES} ${NOT_GMX_SOURCES} ${H5MD_SOURCES})

file(GLOB VOTCA_SOURCES *.cc potentialfunctions/*.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

add_library(votca_csg ${VOTCA_SOURCES} ${IO_SOURCES})
if (TARGET GMX::libgromacs)
  #general to a gmx builds
  target_sources(votca_csg PRIVATE ${GMX_SOURCES})
  target_link_libraries(votca_csg PRIVATE GMX::libgromacs)
endif()


find_package(HDF5 COMPONENTS "C")
if(HDF5_FOUND)
  set(H5MD true)
  target_include_directories (votca_csg PRIVATE ${HDF5_INCLUDE_DIRS})
  target_sources(votca_csg PRIVATE ${H5MD_SOURCES})
  target_link_libraries(votca_csg PRIVATE ${HDF5_LIBRARIES})
else(HDF5_FOUND)
  set(H5MD)
endif(HDF5_FOUND)

configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)

#CMAKE_CURRENT_BINARY_DIR for gitversion.h /votca_config.h
target_include_directories(votca_csg PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_custom_target(gitversion-csg COMMAND ${CMAKE_COMMAND}
  -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)

add_dependencies(votca_csg gitversion-csg)
set_target_properties(votca_csg PROPERTIES SOVERSION ${SOVERSION})
target_link_libraries(votca_csg PUBLIC VOTCA::votca_tools Boost::boost PRIVATE Boost::program_options Boost::filesystem Boost::system Threads::Threads Eigen3::Eigen)
if(MKL_FOUND)
  target_link_libraries(votca_csg PRIVATE MKL::MKL)
endif()

install(TARGETS votca_csg EXPORT VOTCA_CSG_Targets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT VOTCA_CSG_Targets FILE VOTCA_CSG_Targets.cmake NAMESPACE VOTCA:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VOTCA_CSG)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("VOTCA_CSGConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_CSGConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_CSGConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VOTCA_CSG)
add_library(VOTCA::votca_csg ALIAS votca_csg)

