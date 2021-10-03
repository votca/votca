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

set (GIT_HEADER "gitversion.h")
set (NEW_GIT_HEADER "new_gitversion.h")
file(WRITE ${NEW_GIT_HEADER} "#include<string>\nstatic const std::string gitversion = \"${THIS_GIT_ID}\";\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${NEW_GIT_HEADER} ${GIT_HEADER})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${NEW_GIT_HEADER})
