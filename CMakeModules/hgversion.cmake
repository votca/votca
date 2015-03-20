if (GIT_EXECUTABLE)
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
  set(THIS_HG_ID "gitid: ${THIS_GIT_ID}")
elseif (MERCURIAL_EXECUTABLE)
  set(CMAKE_MODULE_PATH ${TOP_SOURCE_DIR}/CMakeModules)
  find_package(Mercurial) #for MERCURIAL_HG_INFO
  MERCURIAL_HG_INFO(${TOP_SOURCE_DIR} THIS)
  MESSAGE("Current merurial revision is ${THIS_HG_ID}")
  set (THIS_HG_ID "hgid: ${THIS_HG_ID}")
else()
  set (THIS_HG_ID)
endif()

set (HG_HEADER "hgversion.h")
set (NEW_HG_HEADER "new_hgversion.h")
file(WRITE ${NEW_HG_HEADER} "static const std::string hgversion = \"${THIS_HG_ID}\";\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${NEW_HG_HEADER} ${HG_HEADER})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${NEW_HG_HEADER})
