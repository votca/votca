cmake_minimum_required(VERSION 2.6)

set(PROJECT_VERSION "1.2.3")

project(csgapps)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)
include(CheckCXXCompilerFlag)

option(EXTERNAL_BOOST "Use external boost" ON)
if (EXTERNAL_BOOST)
  find_package(Boost 1.39.0 COMPONENTS program_options )
  if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost not found,  make sure you have also installed boost and it's dev package or use our internal replacement, if you have built tools with -DEXTERNAL_BOOST=OFF and set it here as well")
  endif(NOT Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  set (BOOST_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
else(EXTERNAL_BOOST)
  find_package(VOTCA_BOOST REQUIRED)
  include_directories(${VOTCA_BOOST_INCLUDE_DIRS})
  set (BOOST_LIBRARIES ${VOTCA_BOOST_LIBRARIES})
endif(EXTERNAL_BOOST)

find_package(VOTCA_TOOLS REQUIRED)
include_directories(${VOTCA_TOOLS_INCLUDE_DIRS})
find_package(VOTCA_CSG REQUIRED)
include_directories(${VOTCA_CSG_INCLUDE_DIRS})

file(GLOB_RECURSE DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} */CMakeLists.txt)
foreach(DIR ${DIRS})
  string(REGEX MATCH ^[^/]* DIR ${DIR})
  file(GLOB SRC ${DIR}/*.cc)
  add_executable(${DIR} ${SRC})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DIR}) #needed for out-of-source build
  target_link_libraries(${DIR} ${VOTCA_CSG_LIBRARIES} ${VOTCA_TOOLS_LIBRARIES} ${BOOST_LIBRARIES})
  install(TARGETS ${DIR} RUNTIME DESTINATION bin)
  set_target_properties(${DIR} PROPERTIES OUTPUT_NAME ${DIR}/csg_${DIR})
endforeach(DIR)
