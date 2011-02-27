file(GLOB_RECURSE BOOST_SOURCES libs/*.cpp)
file(GLOB_RECURSE NOT_BOOST_SOURCES libs/detail/utf8_codecvt_facet.cpp)
list(REMOVE_ITEM BOOST_SOURCES ${NOT_BOOST_SOURCES})

add_library(votca_boost ${BOOST_SOURCES})
install(TARGETS votca_boost LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION ${LIB})

configure_file(libvotca_boost.pc.in libvotca_boost.pc @ONLY)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/libvotca_boost.pc DESTINATION ${LIB}/pkgconfig)

install(DIRECTORY boost DESTINATION include/votca)
