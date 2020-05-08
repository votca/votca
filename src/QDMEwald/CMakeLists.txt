add_subdirectory(src)

if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif(ENABLE_TESTING)
