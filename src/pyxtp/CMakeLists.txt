
find_package(pybind11 REQUIRED)
find_package(VOTCA_XTP REQUIRED NO_MODULE)

file(GLOB_RECURSE BIND_SOURCES binds/*.cc)

pybind11_add_module(xtp_binds ${BIND_SOURCES})
target_include_directories(xtp_binds
 PUBLIC
 binds/include
)

target_compile_options(xtp_binds
 PUBLIC
)
set_target_properties(xtp_binds
  PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
)


target_link_libraries(xtp_binds
 PUBLIC
 VOTCA::votca_xtp
)

# # DPYXTP_VERSION_INFO is defined in the setup.py file
target_compile_definitions(xtp_binds PRIVATE VERSION_INFO=${XTP_BINDS_VERSION_INFO})