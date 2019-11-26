cmake_minimum_required(VERSION 3.10)
cmake_policy(VERSION 3.10)
project(QDMEwald)

option(QDMEwald_USE_CUDA "Using CUDA backend for GPU support" OFF)
option(QDMEwald_KOKKOS_ARCH "Settings for the Kokkos architecture used" "BDW")

option(QDMEwald_EWALD_KONLY "executable only computes k-space contribution of Ewald solver" OFF)
option(QDMEwald_EWALD_RONLY "executable only computes r-space and smaller contributions of Ewald solver" OFF)

set(MAIN_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
#ToDo: figure out how to use these
set(KOKKOS_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/kokkos")
set(CABANA_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/cabana")

if (QDMEwald_USE_CUDA)
    find_package(CUDA)
endif (QDMEwald_USE_CUDA)

set(KOKKOS_ARCH "${QDMEwald_KOKKOS_ARCH}" CACHE STRING "Kokkos architecture string")
set(KOKKOS_ENABLE_SERIAL ON CACHE BOOL "Kokkos serial option")
set(KOKKOS_ENABLE_OPENMP ON CACHE BOOL "Kokkos OpenMP option")

message("${KOKKOS_ARCH} ${KOKKOS_ENABLE_SERIAL} ${KOKKOS_ENABLE_OPENMP}")

include(ExternalProject)

#add OpenMP libraries
find_package(OpenMP REQUIRED)
list(APPEND QDMEwald_LIBS OpenMP::OpenMP_CXX)

#add librt for Kokkos profiling
find_library(LIBRT rt REQUIRED)
list(APPEND QDMEwald_LIBS ${LIBRT})

list(APPEND QDMEwald_LIBS ${CMAKE_DL_LIBS})

#add Kokkos
set(KOKKOS_INSTALL_DIR kokkos_install)

message (
    "KOKKOS passed install dir: ${KOKKOS_INSTALL_PREFIX}"
)
if (QDMEwald_USE_CUDA)
    ExternalProject_add (
        Kokkos_build
        GIT_REPOSITORY https://github.com/kokkos/kokkos.git
        GIT_TAG origin/master
        GIT_SHALLOW TRUE
        SOURCE_DIR ${KOKKOS_INSTALL_DIR}
        #CONFIGURE_COMMAND CXX=${CMAKE_BINARY_DIR}/${KOKKOS_INSTALL_DIR}/bin/nvcc_wrapper cmake <SOURCE_DIR> -DKOKKOS_ARCH=${KOKKOS_ARCH} 
        CONFIGURE_COMMAND cmake <SOURCE_DIR>
                          -DKOKKOS_ARCH=${KOKKOS_ARCH} 
                          -DKOKKOS_ENABLE_SERIAL:BOOL=ON 
                          -DKOKKOS_ENABLE_OPENMP:BOOL=ON
                          -DKOKKOS_ENABLE_CUDA:BOOL=ON
                          -DKOKKOS_ENABLE_CUDA_LAMBDA:BOOL=ON
                          -DKOKKOS_ENABLE_PROFILING:BOOL=OFF
                          -DCMAKE_CXX_COMPILER=${CMAKE_BINARY_DIR}/${KOKKOS_INSTALL_DIR}/bin/nvcc_wrapper
                          -DCMAKE_INSTALL_PREFIX=${KOKKOS_INSTALL_PREFIX}
    )
else (QDMEwald_USE_CUDA)
    ExternalProject_add (
        Kokkos_build
        GIT_REPOSITORY https://github.com/kokkos/kokkos.git
        GIT_TAG origin/master
        GIT_SHALLOW TRUE
        SOURCE_DIR ${KOKKOS_INSTALL_DIR}
        CONFIGURE_COMMAND cmake <SOURCE_DIR> -DKOKKOS_ARCH=${KOKKOS_ARCH} 
                                             -DKOKKOS_ENABLE_SERIAL:BOOL=ON 
                                             -DKOKKOS_ENABLE_OPENMP:BOOL=ON
                                             -DKOKKOS_ENABLE_PROFILING:BOOL=OFF
                                             -DCMAKE_INSTALL_PREFIX=${KOKKOS_INSTALL_PREFIX}
    )
endif (QDMEwald_USE_CUDA)

ExternalProject_get_property (
    Kokkos_build
    INSTALL_DIR
)
message (
    "KOKKOS install dir: ${INSTALL_DIR}"
)
set (
    KOKKOS_BUILD_DIR
    ${INSTALL_DIR}
)
list(APPEND QDMEwald_LIBS ${KOKKOS_INSTALL_PREFIX}/lib/libkokkos.a)
include_directories(${KOKKOS_INSTALL_PREFIX}/include)

#setting Cabana up
#[[
if (QDMEwald_USE_CUDA)
    ExternalProject_add (
        Cabana_build
        GIT_REPOSITORY https://github.com/ECP-copa/Cabana.git
        GIT_TAG origin/master
        GIT_SHALLOW TRUE
        DEPENDS Kokkos_build
        SOURCE_DIR cabana_build
        CONFIGURE_COMMAND cmake <SOURCE_DIR> -DCMAKE_PREFIX_PATH=${KOKKOS_INSTALL_PREFIX}
                                             -DCMAKE_CXX_COMPILER=${CMAKE_BINARY_DIR}/${KOKKOS_INSTALL_DIR}/bin/nvcc_wrapper
                                             -DCabana_ENABLE_Serial:BOOL=ON 
                                             -DCabana_ENABLE_OpenMP:BOOL=ON
                                             -DCabana_ENABLE_Cuda:BOOL=ON
                                             -DCabana_ENABLE_MPI:BOOL=ON
                                             -DCMAKE_INSTALL_PREFIX=${CABANA_INSTALL_PREFIX}
    )
else (QDMEwald_USE_CUDA)
    ExternalProject_add (
        Cabana_build
        GIT_REPOSITORY https://github.com/ECP-copa/Cabana.git
        GIT_TAG origin/master
        GIT_SHALLOW TRUE
        DEPENDS Kokkos_build
        SOURCE_DIR cabana_build
        CONFIGURE_COMMAND cmake <SOURCE_DIR> -DCMAKE_PREFIX_PATH=${KOKKOS_INSTALL_PREFIX}
                                             -DCabana_ENABLE_Serial:BOOL=ON 
                                             -DCabana_ENABLE_OpenMP:BOOL=ON
                                             -DCabana_ENABLE_MPI:BOOL=ON
                                             -DCMAKE_INSTALL_PREFIX=${CABANA_INSTALL_PREFIX}
    )
endif (QDMEwald_USE_CUDA)
list(APPEND QDMEwald_LIBS ${CABANA_INSTALL_PREFIX}/lib64/libcabanacore.a)
include_directories(${CABANA_INSTALL_PREFIX}/include)

#]]

if (QDMEwald_USE_CUDA)
    set(QDMEwald_EXE_NAME "QDMEwald.cuda")
else (QDMEwald_USE_CUDA)
    set(QDMEwald_EXE_NAME "QDMEwald")
endif (QDMEwald_USE_CUDA)


file( GLOB QDMEwald_SRC_FILES src/* )
add_executable(${QDMEwald_EXE_NAME} ${QDMEwald_SRC_FILES})

find_package(MPI)
target_link_libraries(${QDMEwald_EXE_NAME} ${MPI_CXX_LIBRARIES})
include_directories(${MPI_CXX_INCLUDE_DIRS})

set_property(TARGET ${QDMEwald_EXE_NAME} PROPERTY CXX_STANDARD 11)
target_link_libraries(${QDMEwald_EXE_NAME} ${QDMEwald_LIBS})

add_custom_target(
    Kokkos_flags
    VERBATIM
    USES_TERMINAL
    COMMAND grep -i "cflags:" ${KOKKOS_INSTALL_PREFIX}/lib/pkgconfig/kokkos.pc | cut -c 26- | sed "s/;/ /g" > ${CMAKE_BINARY_DIR}/Kokkos.flags
    BYPRODUCTS ${CMAKE_BINARY_DIR}/Kokkos.flags
    DEPENDS Kokkos_build
)

add_dependencies(Kokkos_flags Kokkos_build)
add_dependencies(${QDMEwald_EXE_NAME} Kokkos_build)
add_dependencies(${QDMEwald_EXE_NAME} Kokkos_flags)
#add_dependencies(${QDMEwald_EXE_NAME} Cabana_build)

#execute_process(COMMAND "grep -i \"cflags:\" ${KOKKOS_INSTALL_PREFIX}/lib/pkgconfig/kokkos.pc | cut -c 9-" OUTPUT_VARIABLE KOKKOS_CFLAGS)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} `/bin/cat\ ${CMAKE_BINARY_DIR}/Kokkos.flags`")

message("CXX compiler call: ${CMAKE_CXX_COMPILER}")
message("CXX CFLAGS: ${CMAKE_CXX_FLAGS}")
message("Kokkos CFLAGS: ${KOKKOS_CFLAGS}")

# collect compile flags for the different CMake options
if (QDMEwald_USE_CUDA)
    set( QDMEwald_CXXFLAGS "${QDMEwald_CXXFLAGS} --expt-extended-lambda -w -DCMAKE_KOKKOS_CUDA_ENABLED") 
endif (QDMEwald_USE_CUDA)
if (QDMEwald_EWALD_KONLY)
    set( QDMEwald_CXXFLAGS "${QDMEwald_CXXFLAGS} -DCMAKE_QDMEwald_EWALD_KONLY") 
endif (QDMEwald_EWALD_KONLY)
if (QDMEwald_EWALD_RONLY)
    set( QDMEwald_CXXFLAGS "${QDMEwald_CXXFLAGS} -DCMAKE_QDMEwald_EWALD_RONLY") 
endif (QDMEwald_EWALD_RONLY)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${QDMEwald_CXXFLAGS}" )

if(ENABLE_TESTING)
  enable_testing()
  include(FindPackageHandleStandardArgs) 
  add_subdirectory(tests)
endif(ENABLE_TESTING)
