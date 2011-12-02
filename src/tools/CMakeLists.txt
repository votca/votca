add_executable(ctp_map ctp_map.cc md2qm_observer.cc)
target_link_libraries(ctp_map votca_ctp ${BOOST_LIBRARIES})
add_executable(ctp_run ctp_run.cc)
target_link_libraries(ctp_run votca_ctp ${BOOST_LIBRARIES})
add_executable(ctp_test ctp_test.cc)
target_link_libraries(ctp_test votca_ctp ${BOOST_LIBRARIES})

install(TARGETS ctp_map ctp_run ctp_test RUNTIME DESTINATION bin)
