add_custom_command(OUTPUT ${CURRENT_DIR}/csgid.tex
    COMMAND VOTCASHARE=${CSG_SHARE} $<TARGET_FILE:VOTCA::csg_call> --help > csg_call.help
    COMMAND sed -n -e "'s/.*version \\(.*\\)/\\1/p'" csg_call.help > csg_call.version
    COMMAND sed -e "'s/_/\\\\&/g'" -e "'s/.*/\\\\newcommand{\\\\csgid}{&}/'" csg_call.version > ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/csgid.tex
    DEPENDS VOTCA::csg_call)
