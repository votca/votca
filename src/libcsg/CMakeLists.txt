file(GLOB GMX_SOURCES modules/io/gmx*.cc)

file(GLOB H5MD_SOURCES modules/io/h5md*.cc)

file(GLOB IO_SOURCES modules/io/*.cc)
list(REMOVE_ITEM IO_SOURCES ${GMX_SOURCES} ${NOT_GMX_SOURCES} ${H5MD_SOURCES})

file(GLOB VOTCA_SOURCES *.cc potentialfunctions/*.cc)
file(GLOB NOT_VOTCA_SOURCES version_nb.cc test.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

add_library(votca_csg ${VOTCA_SOURCES} ${IO_SOURCES})
if (WITH_GMX)
  #general to a gmx builds
  target_sources(votca_csg PRIVATE ${GMX_SOURCES})
  target_link_libraries(votca_csg PRIVATE GMX::libgromacs)
endif(WITH_GMX)


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
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(gitversion-csg COMMAND ${CMAKE_COMMAND}
  -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
  -P ${CMAKE_MODULE_PATH}/gitversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES gitversion.h)

add_dependencies(votca_csg gitversion-csg)
set_target_properties(votca_csg PROPERTIES SOVERSION ${SOVERSION})
target_link_libraries(votca_csg PUBLIC ${VOTCA_TOOLS_LIBRARIES} ${Boost_LIBRARIES} PRIVATE ${MKL_LIBRARIES} ${THREAD_LIBRARIES})
install(TARGETS votca_csg LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
