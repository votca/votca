cmake_minimum_required(VERSION 3.12)

project(pyxtp)

find_package(pybind11 REQUIRED)
find_package(VOTCA_XTP REQUIRED NO_MODULE)

pybind11_add_module(pyxtp binds/pyxtp.cc)


target_link_libraries(pyxtp PUBLIC VOTCA::votca_xtp)

# DPYXTP_VERSION_INFO is defined in the setup.py file
target_compile_definitions(pyxtp PRIVATE VERSION_INFO=${DPYXTP_VERSION_INFO})