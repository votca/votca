enable_language(C)
include(CheckCCompilerFlag)

add_library(votca_expat libs/xmlparse.c libs/xmlrole.c libs/xmltok.c)
install(TARGETS votca_expat LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

add_definitions( -DHAVE_EXPAT_CONFIG_H )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libvotca_expat.pc.in 
  ${CMAKE_CURRENT_BINARY_DIR}/libvotca_expat.pc)
