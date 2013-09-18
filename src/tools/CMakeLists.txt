add_executable(ctp_map ctp_map.cc Md2QmEngine.cc)
target_link_libraries(ctp_map votca_ctp ${BOOST_LIBRARIES})
add_executable(ctp_run ctp_run.cc)
target_link_libraries(ctp_run votca_ctp ${BOOST_LIBRARIES})
add_executable(ctp_app ctp_app.cc)
target_link_libraries(ctp_app votca_ctp ${BOOST_LIBRARIES})

install(TARGETS ctp_map ctp_run ctp_app RUNTIME DESTINATION bin)
