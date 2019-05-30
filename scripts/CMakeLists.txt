option(WITH_RC_FILES "Install votca rc files, no need when installing under /usr" ON)
if (WITH_RC_FILES)
  if(APPLE)
    set(LD_LIBRARY_PATH "DYLD_LIBRARY_PATH")
  else(APPLE)
    include(CheckCSourceCompiles)
    check_c_source_compiles(
      "int main(void) {
      #ifndef _AIX
      #error
      #endif
      return 0;
      }" HAVE_AIX)
    if(HAVE_AIX)
      set(LD_LIBRARY_PATH "LIBPATH")
    else(HAVE_AIX)
      set(LD_LIBRARY_PATH "LD_LIBRARY_PATH")
    endif(HAVE_AIX)
  endif(APPLE)
  configure_file(VOTCARC.csh.in ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.csh @ONLY)
  configure_file(VOTCARC.bash.in ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.bash @ONLY)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.csh
    ${CMAKE_CURRENT_BINARY_DIR}/VOTCARC.bash DESTINATION ${CMAKE_INSTALL_BINDIR})
endif (WITH_RC_FILES)

foreach(SCRIPT votca_compare)
  configure_file(${SCRIPT}.in ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.out @ONLY)
  add_custom_target(${SCRIPT}_build ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT})
  add_custom_command(OUTPUT ${SCRIPT} COMMAND ${CMAKE_COMMAND}
    -DINPUT="${SCRIPT}.out" -DOUTPUT="${SCRIPT}"
    -DGIT_EXECUTABLE="${GIT_EXECUTABLE}"
    -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${CMAKE_MODULE_PATH}/gitscript.cmake)
  if(ENABLE_TESTING)
    add_custom_command(TARGET ${SCRIPT}_build POST_BUILD COMMAND chmod +x ${SCRIPT})
  endif(ENABLE_TESTING)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${SCRIPT})
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach(SCRIPT)

add_executable(VOTCA::votca_compare IMPORTED GLOBAL)
set_property(TARGET VOTCA::votca_compare PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/votca_compare")
