add_custom_command(OUTPUT ${CURRENT_DIR}/config.t2t
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/../config.t2t ${CURRENT_DIR}/config.t2t
  DEPENDS ${CURRENT_DIR}/../config.t2t)

add_custom_command(OUTPUT ${CURRENT_DIR}/csg_table.tex
  COMMAND VOTCASHARE=${CSG_SHARE} VOTCA_CSG_DEFAULTS=${VOTCA_CSG_DEFAULTS} ${CSG_CALL} -l > ${CURRENT_DIR}/csg_table.txt
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/call2t2t.sh ${CURRENT_DIR}/csg_table.txt
    > ${CURRENT_DIR}/csg_table.t2t
  COMMAND ${TXT2TAGS_EXECUTABLE} -t tex --no-header --infile ${CURRENT_DIR}/csg_table.t2t
      --outfile ${CURRENT_DIR}/csg_table.tex
      DEPENDS ${CSG_CALL} ${CURRENT_DIR}/config.t2t)

#don't call "csg_call -l" as csg_call isn't ready at cmake time in case of a super-build
file(READ ${CSG_SHARE}/scripts/inverse/csg_table CSG_SCRIPTS)
string(REGEX REPLACE "#[^\n]*\n" "" CSG_SCRIPTS "${CSG_SCRIPTS}") #remove comments
string(REGEX REPLACE "\n" ";" CSG_SCRIPTS "${CSG_SCRIPTS}") #turn into list
string(REGEX REPLACE "; *;" ";" CSG_SCRIPTS "${CSG_SCRIPTS}") #remove empty line
string(REGEX REPLACE "[^ ;]+ +[^ ;]+ +([^ ;]+) *[^ ;]*[^;]*;" "\\1;" CSG_SCRIPTS "${CSG_SCRIPTS}")
list(SORT CSG_SCRIPTS)
list(REMOVE_DUPLICATES CSG_SCRIPTS)
set(SCRIPT_TEXS)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake)
foreach(script ${CSG_SCRIPTS}) 
  string(TOUPPER "${script}" SCRIPT)
  add_custom_command(OUTPUT ${CURRENT_DIR}/${script}.tex
    COMMAND VOTCASHARE=${CSG_SHARE} VOTCA_CSG_DEFAULTS=${VOTCA_CSG_DEFAULTS} CSG_CALL=${CSG_CALL} ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/script2t2t.sh ${script} > ${CURRENT_DIR}/${script}.t2t
    COMMAND ${TXT2TAGS_EXECUTABLE} -t tex --no-header --infile ${CURRENT_DIR}/${script}.t2t
      --outfile ${CURRENT_DIR}/${script}.tex
    DEPENDS ${CSG_SHARE}/scripts/inverse/${script} ${CURRENT_DIR}/config.t2t)
  list(APPEND SCRIPT_TEXS ${CURRENT_DIR}/${script}.tex)
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake "\\input{${CURRENT_DIR}/${script}}\n")
endforeach()
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex DEPENDS ${PROG_TEXS}
  DEPENDS ${SCRIPT_TEXS})
