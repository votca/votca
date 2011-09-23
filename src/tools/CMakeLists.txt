add_executable(ctp_map ctp_map.cc md2qm_observer.cc)
target_link_libraries(ctp_map votca_ctp)
add_executable(ctp_run ctp_run.cc)
target_link_libraries(ctp_run votca_ctp)
add_executable(ctp_polymer polymer.cc)
target_link_libraries(ctp_polymer votca_ctp)

install(TARGETS ctp_map ctp_run ctp_polymer RUNTIME DESTINATION bin)
