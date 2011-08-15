cmake_minimum_required(VERSION 2.6)

set(PROJECT_VERSION "1.3.0-dev")

project(csgapps)

file(GLOB_RECURSE DIRS */CMakeLists.txt)
foreach(DIR ${DIRS})
  string(REGEX MATCH ^.*/ DIR ${DIR})
  add_subdirectory(${DIR})
  unset(NAME CACHE)
endforeach(DIR)
