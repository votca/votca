if (GIT_EXECUTABLE)
  # refresh git index
  execute_process(COMMAND ${GIT_EXECUTABLE} update-index -q --refresh
    WORKING_DIRECTORY ${TOP_SOURCE_DIR} TIMEOUT 5 OUTPUT_QUIET
    ERROR_VARIABLE EXEC_ERR OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  #later use git describe here
  execute_process( COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${TOP_SOURCE_DIR}
    OUTPUT_VARIABLE THIS_GIT_ID OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process( COMMAND ${GIT_EXECUTABLE} diff-index --name-only HEAD
    WORKING_DIRECTORY ${TOP_SOURCE_DIR}
    OUTPUT_VARIABLE _HAS_CHANGES OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if (NOT "${_HAS_CHANGES}" STREQUAL "")
    set(THIS_GIT_ID "${THIS_GIT_ID} (dirty)")
  endif()
  message("Current git revision is ${THIS_GIT_ID}")
  set(THIS_GIT_ID "gitid: ${THIS_GIT_ID}")
else()
  set (THIS_GIT_ID)
endif()
file(READ ${INPUT} CONTENT)
string(REGEX REPLACE "#CSG_GIT_ID#" "${THIS_GIT_ID}" NEW_CONTENT "${CONTENT}")
file(WRITE "${OUTPUT}.tmp" "${NEW_CONTENT}")
file(COPY "${OUTPUT}.tmp" DESTINATION "${OUTPUT}.tmp.x" FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${OUTPUT}.tmp.x/${OUTPUT}.tmp ${OUTPUT})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${OUTPUT}.tmp ${OUTPUT}.tmp.x/${OUTPUT}.tmp)
execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${OUTPUT}.tmp.x)
