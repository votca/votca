cmake_minimum_required(VERSION 3.12)

project(pyxtp)

find_package(pybind11 REQUIRED)
find_package(VOTCA_XTP REQUIRED NO_MODULE)

pybind11_add_module(xtp_binds binds/pyxtp.cc)


target_link_libraries(xtp_binds PUBLIC VOTCA::votca_xtp)

# DPYXTP_VERSION_INFO is defined in the setup.py file
target_compile_definitions(xtp_binds PRIVATE VERSION_INFO=${DXTP_BINDS_VERSION_INFO})