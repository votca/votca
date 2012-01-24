add_executable(ctp_map ctp_map.cc md2qm_observer.cc)
target_link_libraries(ctp_map votca_ctp ${BOOST_LIBRARIES})
add_executable(ctp_run ctp_run.cc)
target_link_libraries(ctp_run votca_ctp ${BOOST_LIBRARIES})
add_executable(ctp_map_exp ctp_map_exp.cc Md2QmEngine.cc)
target_link_libraries(ctp_map_exp votca_ctp ${BOOST_LIBRARIES})

install(TARGETS ctp_map ctp_run ctp_map_exp RUNTIME DESTINATION bin)
