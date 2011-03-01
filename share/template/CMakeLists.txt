add_custom_command(OUTPUT votca/csg
 COMMAND ${CMAKE_COMMAND} -E make_directory votca
 COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/include votca/csg
)
add_custom_target(csg_include_links DEPENDS votca/csg)
set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM TRUE)

foreach(PROG template template_threaded)
  file(GLOB ${PROG}_SOURCES ${PROG}.cc)
  add_executable(${PROG} ${${PROG}_SOURCES})
  add_dependencies(${PROG} csg_include_links)
  target_link_libraries(${PROG} votca_csg)
 install(FILES ${${PROG}_SOURCES} DESTINATION @DATA@/template)
endforeach(PROG)

install(FILES Makefile.template DESTINATION @DATA@/template RENAME Makefile)
install(FILES CMakeLists.txt.template DESTINATION @DATA@/template RENAME CMakeLists.txt)
