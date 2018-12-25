add_custom_command(OUTPUT ${CURRENT_DIR}/csgid.tex
    COMMAND VOTCASHARE=${CSG_SHARE} ${CSG_CALL} --help > csg_call.help
    COMMAND sed -n -e "'s/.*version \\(.*\\)/\\\\newcommand{\\\\csgid}{\\1}/p'" csg_call.help
      > ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/csgid.tex
    DEPENDS ${CSG_CALL})
