# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)


########################################################################
# Checks what linear algebra packages are installed                    #
########################################################################
option(WITH_LINALG "Checking which BLAS libraries are installed" ON)
if (WITH_LINALG)

  find_package(Eigen)
  find_package(GSL)
  find_package(MKL)

  if (EIGEN_FOUND)
    message("-- EIGEN found (not used at the moment)")
    set(WITH_EIGEN TRUE) #used in votca_ctp_config.h
    include_directories(${EIGEN_INCLUDE_DIRS})
  else()
    message("-- EIGEN not found")
    set(WITH_EIGEN FALSE) #used in votca_ctp_config.h
  endif()

  if (GSL_FOUND)
    message("-- GSL found")
    include_directories(${GSL_INCLUDE_DIRS})
    set(GSL_PKG "gsl")
    set(WITH_GSL TRUE) #used in votca_ctp_config.h
  else()
    message(WARNING "-- GSL not found")
    set(WITH_GSL FALSE) #used in votca_ctp_config.h
  endif()

  if (MKL_FOUND)
    message("-- MKL found")
    set(WITH_MKL TRUE) #used in votca_ctp_config.h
  else()
    message("-- MKL not found")
    set(WITH_MKL FALSE) #used in votca_ctp_config.h
  endif()

  if (WITH_MKL AND CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message("-- MKL wrappers are used for linear algebra")
    file(GLOB VOTCA_LINALG_MKL_SOURCES mkl/*.cc)
    set(VOTCA_LINALG_SOURCES ${VOTCA_LINALG_MKL_SOURCES} PARENT_SCOPE)
  else (WITH_GSL)
    message("-- GSL wrappers are used for linear algebra")
    file(GLOB  VOTCA_LINALG_GSL_SOURCES gsl/*.cc)
    set(VOTCA_LINALG_SOURCES ${VOTCA_LINALG_GSL_SOURCES} PARENT_SCOPE)
  endif()

endif(WITH_LINALG)
