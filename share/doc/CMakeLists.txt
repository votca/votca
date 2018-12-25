set(DOXYGEN_SOURCES)
foreach(DIR ${ENABLED_VOTCA_PACKAGES})
  set(ABS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../${DIR})
  set(INPUT_DIRS "${INPUT_DIRS} ${ABS_DIR}")
  file(GLOB_RECURSE ${DIR}_SOURCES ${ABS_DIR}/*.cc ${ABS_DIR}/*.h
         ${ABS_DIR}/*.dox ${ABS_DIR}/guide/*.md)
  list(APPEND DOXYGEN_SOURCES ${${DIR}_SOURCES})
endforeach()
file(GLOB EXTRA_DOXYGEN_PAGES ${CMAKE_CURRENT_SOURCE_DIR}/*.md ${CMAKE_SOURCE_DIR}/*.md)
list(APPEND DOXYGEN_SOURCES ${EXTRA_DOXYGEN_PAGES})
string(REPLACE ";" " " INPUT_SOURCES "${DOXYGEN_SOURCES}")
configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doxygen DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
  COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile ${DOXYGEN_SOURCES}
  COMMENT "Build doxygen documentation")
