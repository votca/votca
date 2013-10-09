add_executable(ctp_map ctp_map.cc Md2QmEngine.cc)
target_link_libraries(ctp_map votca_ctp ${BOOST_LIBRARIES})

add_executable(ctp_run ctp_run.cc)
target_link_libraries(ctp_run votca_ctp ${BOOST_LIBRARIES})

add_executable(ctp_tools ctp_tools.cc)
target_link_libraries(ctp_tools votca_ctp ${BOOST_LIBRARIES})

add_executable(ctp_evaluate ctp_evaluate.cc)
target_link_libraries(ctp_evaluate votca_ctp ${BOOST_LIBRARIES})

add_executable(ctp_dump ctp_dump.cc)
target_link_libraries(ctp_dump votca_ctp ${BOOST_LIBRARIES})

install(TARGETS ctp_map ctp_tools ctp_evaluate ctp_run ctp_dump RUNTIME DESTINATION bin)
