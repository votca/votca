add_custom_command(OUTPUT ${CURRENT_DIR}/config.t2t
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/../config.t2t ${CURRENT_DIR}/config.t2t
  DEPENDS ${CURRENT_DIR}/../config.t2t)

set(PROG_TEXS)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake)
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
  add_custom_command(OUTPUT ${CURRENT_DIR}/${prog}.tex
    COMMAND ${${PROG}} --help > ${CURRENT_DIR}/${prog}.help
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/help2t2t.sh ${prog} ${CURRENT_DIR}/${prog}.help
      > ${CURRENT_DIR}/${prog}.t2t
    COMMAND ${TXT2TAGS_EXECUTABLE} -t tex --no-header --infile ${CURRENT_DIR}/${prog}.t2t
      --outfile ${CURRENT_DIR}/${prog}.tex
      DEPENDS ${${PROG}} ${CURRENT_DIR}/config.t2t)
  list(APPEND PROG_TEXS ${CURRENT_DIR}/${prog}.tex)
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake "\\input{${CURRENT_DIR}/${prog}}\n")
endforeach()

set(VOTCASHARE ${CMAKE_CURRENT_SOURCE_DIR}/../share)
foreach(PROG ${XTP_SCRIPTS})
  get_filename_component(prog ${PROG} NAME_WE)
  add_custom_command(OUTPUT ${CURRENT_DIR}/${prog}.tex
    COMMAND chmod +x ${PROG}
    COMMAND VOTCASHARE=${VOTCASHARE} ${PROG} --help > ${CURRENT_DIR}/${prog}.help
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/help2t2t.sh ${prog} ${CURRENT_DIR}/${prog}.help
      > ${CURRENT_DIR}/${prog}.t2t
    COMMAND ${TXT2TAGS_EXECUTABLE} -t tex --no-header --infile ${CURRENT_DIR}/${prog}.t2t
      --outfile ${CURRENT_DIR}/${prog}.tex
      DEPENDS ${PROG} ${CURRENT_DIR}/config.t2t)
  list(APPEND PROG_TEXS ${CURRENT_DIR}/${prog}.tex)
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake "\\input{${CURRENT_DIR}/${prog}}\n")
endforeach()

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex 
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex DEPENDS ${PROG_TEXS})
