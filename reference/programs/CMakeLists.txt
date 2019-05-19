add_custom_command(OUTPUT ${CURRENT_DIR}/config.t2t
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/../config.t2t ${CURRENT_DIR}/config.t2t
  DEPENDS ${CURRENT_DIR}/../config.t2t)

set(PROG_TEXS)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake)
foreach(prog csg_boltzmann csg_call csg_density csg_dlptopol csg_dump csg_fmatch csg_gmxtopol csg_imcrepack csg_inverse csg_map
  csg_property csg_resample csg_reupdate csg_stat)
  add_custom_command(OUTPUT ${CURRENT_DIR}/${prog}.tex
    COMMAND VOTCASHARE=${CSG_SHARE} $<TARGET_FILE:VOTCA::${prog}> --help > ${CURRENT_DIR}/${prog}.help
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_DIR}/help2t2t.sh ${prog} ${CURRENT_DIR}/${prog}.help
      > ${CURRENT_DIR}/${prog}.t2t
    COMMAND ${TXT2TAGS_EXECUTABLE} -t tex --no-header --infile ${CURRENT_DIR}/${prog}.t2t
      --outfile ${CURRENT_DIR}/${prog}.tex
      DEPENDS ${prog_dep} ${CURRENT_DIR}/config.t2t)
  list(APPEND PROG_TEXS ${CURRENT_DIR}/${prog}.tex)
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake "\\input{${CURRENT_DIR}/${prog}}\n")
endforeach()
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex 
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_DIR}/all.tex DEPENDS ${PROG_TEXS})
