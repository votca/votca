option(WITH_VOTCA_KMCOLD "Enable the calculator which uses the old kmc code!" OFF)
#also used in votca_config.h
if (WITH_VOTCA_KMCOLD)
  find_package(VOTCA_KMCOLD)
  if(NOT VOTCA_KMCOLD_FOUND)
    message(FATAL_ERROR "libvotca_kmcold not found (obsolate package). Please install or better use the new kmc code.")
  endif(NOT VOTCA_KMCOLD_FOUND)
  include_directories(${VOTCA_KMCOLD_INCLUDE_DIRS})
  set(VOTCA_KMCOLD_PKG "libvotca_kmcold")
else(WITH_VOTCA_KMCOLD)
  set(VOTCA_KMCOLD_PKG)
endif(WITH_VOTCA_KMCOLD)


configure_file(votca_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/votca_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(hgversion COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${CMAKE_MODULE_PATH}/hgversion.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES hgversion.h)


file(GLOB_RECURSE VOTCA_SOURCES *.cc)
file(GLOB_RECURSE NOT_VOTCA_SOURCES version_nb.cc)
list(REMOVE_ITEM VOTCA_SOURCES ${NOT_VOTCA_SOURCES})

add_library(votca_ctp  ${VOTCA_SOURCES})
add_dependencies(votca_ctp hgversion)
target_link_libraries(votca_ctp ${VOTCA_CSG_LIBRARIES} ${VOTCA_KMCOLD_LIBRARIES} ${VOTCA_MOO_LIBRARIES} ${VOTCA_TOOLS_LIBRARIES})
install(TARGETS votca_ctp LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_ctp.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libvotca_ctp.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libvotca_ctp.pc DESTINATION ${LIB}/pkgconfig)
