add_executable(ctp_map ctp_map.cc Md2QmEngine.cc)
target_link_libraries(ctp_map votca_ctp ${BOOST_LIBRARIES})

add_executable(ctp_run ctp_run.cc)
target_link_libraries(ctp_run votca_ctp ${BOOST_LIBRARIES})

add_executable(ctp_tools ctp_tools.cc)
target_link_libraries(ctp_tools votca_ctp ${BOOST_LIBRARIES})

add_executable(ctp_parallel ctp_parallel.cc)
target_link_libraries(ctp_parallel votca_ctp ${BOOST_LIBRARIES})

add_executable(ctp_dump ctp_dump.cc)
target_link_libraries(ctp_dump votca_ctp ${BOOST_LIBRARIES})

install(TARGETS ctp_map ctp_dump ctp_tools ctp_run ctp_parallel RUNTIME DESTINATION bin)
