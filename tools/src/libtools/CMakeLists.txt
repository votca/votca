find_package(EXPAT REQUIRED)
set_package_properties(EXPAT PROPERTIES TYPE REQUIRED PURPOSE "Used for XML parsing")

file(GLOB VOTCA_SOURCES *.cc)

add_library(votca_tools ${VOTCA_SOURCES} ${VOTCA_LINALG_SOURCES})
target_compile_features(votca_tools PUBLIC cxx_std_17)

# CMAKE_CURRENT_BINARY_DIR for gitversion.h
# and ${PROJECT_BINARY_DIR}/include/votca/tools for votca_tools_config.h
# and ${PROJECT_BINARY_DIR}/include for votca_tools_config.h when used by
# csg/xtp
target_include_directories(votca_tools PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tools/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/tools/include/votca/tools>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/tools/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_custom_target(gitversion COMMAND ${CMAKE_COMMAND}
  -DTOP_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${PROJECT_SOURCE_DIR}/CMakeModules/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)

add_dependencies(votca_tools gitversion)

set_target_properties(votca_tools PROPERTIES SOVERSION ${SOVERSION})
target_link_libraries(votca_tools PUBLIC Boost::boost Boost::program_options Eigen3::Eigen PRIVATE
	Boost::system EXPAT::EXPAT Threads::Threads ${MATH_LIBRARIES})
if(MKL_FOUND)
  target_link_libraries(votca_tools PUBLIC MKL::MKL)
endif()

if(OpenMP_FOUND)
  target_link_libraries(votca_tools PUBLIC OpenMP::OpenMP_CXX)
endif()

if(FFTW3_FOUND)
  # fftw3.h gets included in our public eigen.h through Eigen/FFT
  target_link_libraries(votca_tools PUBLIC FFTW3::fftw3)
  install(FILES "${PROJECT_SOURCE_DIR}/CMakeModules/FindFFTW3.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VOTCA_TOOLS)
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
