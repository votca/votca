if(WITH_DIPRO)
  enable_language(Fortran)
  find_package(LAPACK REQUIRED) 
  
  add_executable(mergeorbitals.x mergeorbitals.f90 stringmod.f90 ranking.f90)
  add_executable(DIPRO.x DIPRO-TI.f90 stringmod.f90 output.f90 ranking.f90 precmod.f90 parsetools.f90 parsetools_noSCF.f90)
  target_link_libraries(DIPRO.x ${LAPACK_LIBRARIES})
  
  install(TARGETS mergeorbitals.x DIPRO.x RUNTIME DESTINATION bin)
endif(WITH_DIPRO)
