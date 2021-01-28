cmake_minimum_required(VERSION 3.12)

project(pyxtp LANGUAGES CXX)

find_package(pybind11 REQUIRED HINTS)
find_package(VOTCA_XTP REQUIRED NO_MODULE)

# set(MKL_THREAD_LAYER "Intel OpenMP" CACHE STRING "The thread layer to choose for MKL")
# find_package(MKL)
# # https://software.intel.com/en-us/articles/symbol-lookup-error-when-linking-intel-mkl-with-gcc-on-ubuntu
# if(MKL_FOUND AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
# endif()

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


pybind11_add_module(xtp_binds binds/pyxtp.cc)
target_include_directories(xtp_binds
 PUBLIC
 /usr/include
 ${PROJECT_SOURCE_DIR}/binds/include
)

target_compile_options(xtp_binds
 PUBLIC
 "-Wall"
)
set_target_properties(xtp_binds
  PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
)


target_link_libraries(xtp_binds
 PUBLIC
 VOTCA::votca_xtp
 "-Wl,--no-as-needed"
)

# add_executable(main binds/main.cc binds/pyxtp.cc)
# target_include_directories(main
#  PUBLIC
#  ${PROJECT_SOURCE_DIR}/binds/include
#  )

# target_link_libraries(main
#  PUBLIC
#  VOTCA::votca_xtp
#  PRIVATE
#  pybind11::embed
#  )

# # DPYXTP_VERSION_INFO is defined in the setup.py file
target_compile_definitions(xtp_binds PRIVATE VERSION_INFO=${XTP_BINDS_VERSION_INFO})