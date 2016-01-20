foreach(PROG xtp_map xtp_run xtp_tools xtp_parallel xtp_dump)
  file(GLOB ${PROG}_SOURCES ${PROG}*.cc)

  if( ${PROG} STREQUAL "xtp_map" )
     add_executable(xtp_map xtp_map.cc Md2QmEngine.cc)
  else()
     add_executable(${PROG} ${${PROG}_SOURCES})
  endif()

  target_link_libraries(${PROG} votca_xtp ${BOOST_LIBRARIES})
  install(TARGETS ${PROG} RUNTIME DESTINATION bin)

  if (BUILD_MANPAGES)
     add_custom_command(OUTPUT ${PROG}.man
       COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROG} --man > ${PROG}.man
       COMMENT "Building ${PROG} manpage"
       DEPENDS ${PROG})
     add_custom_target(${PROG}_manpage DEPENDS ${PROG}.man)
     add_dependencies(manpages ${PROG}_manpage)
     install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROG}.man DESTINATION ${MAN}/man1 RENAME ${PROG}.1)
  endif (BUILD_MANPAGES)
endforeach(PROG)

