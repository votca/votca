foreach(SCRIPT csg.rc.bash csg.rc.csh  csg_call  csg_inverse help2t2t multi_g_rdf)
  configure_file(${SCRIPT}.in ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT} @ONLY)
endforeach(SCRIPT)
