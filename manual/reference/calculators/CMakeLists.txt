set(XTP_PROGS_FULL)
foreach(prog ${XTP_PROGS})
  string(TOUPPER "${prog}" PROG)
  if(TARGET ${prog})
    set(${PROG} $<TARGET_FILE:${prog}>)
  endif()
  find_program(${PROG} ${prog})
  find_package_handle_standard_args(${PROG} REQUIRED_VARS ${PROG})
  if(NOT ${PROG})
    message(FATAL_ERROR "Could not find ${prog}")
  endif()
  set(XTP_PROGS_FULL "${XTP_PROGS_FULL} ${${PROG}}")
endforeach()

set(XTP_SHARE ${CMAKE_CURRENT_SOURCE_DIR}/../share)
add_custom_command(OUTPUT ${CURRENT_DIR}/calculators.tex
  COMMAND XTP_SHARE=${XTP_SHARE} ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/extract_calculators.sh ${XTP_PROGS_FULL}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}
  DEPENDS ${XTP_PROGS_FULL})
