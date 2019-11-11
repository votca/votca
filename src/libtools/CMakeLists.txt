find_package(EXPAT REQUIRED)

file(GLOB VOTCA_SOURCES *.cc)

add_library(votca_tools ${VOTCA_SOURCES} ${VOTCA_LINALG_SOURCES})

#CMAKE_CURRENT_BINARY_DIR for gitversion.h
#PROJECT_BINARY_DIR/include for votca_config.h
target_include_directories(votca_tools PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_custom_target(gitversion COMMAND ${CMAKE_COMMAND}
  -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)

add_dependencies(votca_tools gitversion)

set_target_properties(votca_tools PROPERTIES SOVERSION ${SOVERSION})
target_link_libraries(votca_tools PUBLIC Boost::boost Boost::program_options Eigen3::Eigen PRIVATE
	Boost::filesystem Boost::system EXPAT::EXPAT Threads::Threads ${MATH_LIBRARIES})
if(TARGET MKL::MKL)
  target_link_libraries(votca_tools PUBLIC MKL::MKL)
  foreach(_PROP IMPORTED_LOCATION INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES IMPORTED_NO_SONAME)
    get_target_property(MKL_${_PROP} MKL::MKL ${_PROP})
  endforeach()
  foreach(_TGT Interface ThreadLayer)
    get_target_property(MKL_${_TGT}_IMPORTED_LOCATION MKL::${_TGT} IMPORTED_LOCATION)
  endforeach()
endif()
if(TARGET FFTW3::FFTW3)
  target_link_libraries(votca_tools PRIVATE FFTW3::FFTW3)
endif()
install(TARGETS votca_tools EXPORT VOTCA_TOOLS_Targets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT VOTCA_TOOLS_Targets FILE VOTCA_TOOLS_Targets.cmake NAMESPACE VOTCA:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VOTCA_TOOLS)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("VOTCA_TOOLSConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
configure_file(VOTCA_TOOLSConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/VOTCA_TOOLSConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_TOOLSConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/VOTCA_TOOLSConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VOTCA_TOOLS)
add_library(VOTCA::votca_tools ALIAS votca_tools)
