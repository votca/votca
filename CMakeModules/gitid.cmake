set(CMAKE_MODULE_PATH ${TOP_SOURCE_DIR}/CMakeModules)

if (GIT_EXECUTABLE)
  # refresh git index
  execute_process(COMMAND ${GIT_EXECUTABLE} update-index -q --refresh
    WORKING_DIRECTORY ${TOP_SOURCE_DIR} TIMEOUT 5 OUTPUT_QUIET
    ERROR_VARIABLE EXEC_ERR OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  #later use git describe here
  execute_process( COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${TOP_SOURCE_DIR}
    OUTPUT_VARIABLE MANUAL_GIT_ID OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process( COMMAND ${GIT_EXECUTABLE} diff-index --name-only HEAD
    WORKING_DIRECTORY ${TOP_SOURCE_DIR}
    OUTPUT_VARIABLE _HAS_CHANGES OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if (NOT "${_HAS_CHANGES}" STREQUAL "")
    set(MANUAL_GIT_ID "${MANUAL_GIT_ID} (dirty)")
  endif()
  message("Current git revision is ${MANUAL_GIT_ID}")
  set(MANUAL_GIT_ID "gitid: ${MANUAL_GIT_ID}")
  set(MANUAL_GIT_ID "(MANUAL_GIT_ID)")
else()
  set (MANUAL_GIT_ID)
endif()

set (GIT_HEADER "gitid.tex")
set (NEW_GIT_HEADER "new_gitid.tex")
string(REPLACE "_" "\\_" SANITIZED_PROJECT_VERSION "${PROJECT_VERSION}")
file(WRITE ${NEW_GIT_HEADER} "\\newcommand{\\gitid}{${SANITIZED_PROJECT_VERSION} ${MANUAL_GIT_ID}}\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${NEW_GIT_HEADER} ${GIT_HEADER})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${NEW_GIT_HEADER})
