# Apple clang ships no OpenMP runtime; Homebrew's libomp is keg-only, so
# FindOpenMP finds it only with OpenMP_ROOT. Set OpenMP_ROOT from
# `brew --prefix libomp` unless the user points elsewhere.
#
# Include this before the first find_package(OpenMP) -- also before packages
# whose config files look for OpenMP themselves (Homebrew's GROMACS does, via
# find_dependency): a failed OpenMP search is cached, so setting OpenMP_ROOT
# afterwards comes too late.
if(APPLE AND NOT DEFINED OpenMP_ROOT AND NOT DEFINED ENV{OpenMP_ROOT})
  find_program(BREW_EXECUTABLE brew)
  mark_as_advanced(BREW_EXECUTABLE)
  if(BREW_EXECUTABLE)
    execute_process(COMMAND ${BREW_EXECUTABLE} --prefix libomp
      OUTPUT_VARIABLE _votca_libomp_prefix OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _votca_libomp_result ERROR_QUIET)
    if(_votca_libomp_result EQUAL 0 AND IS_DIRECTORY "${_votca_libomp_prefix}")
      set(OpenMP_ROOT "${_votca_libomp_prefix}")
      message(STATUS "Using Homebrew libomp: ${OpenMP_ROOT}")
    endif()
    unset(_votca_libomp_prefix)
    unset(_votca_libomp_result)
  endif()
endif()
