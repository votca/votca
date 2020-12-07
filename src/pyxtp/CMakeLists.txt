cmake_minimum_required(VERSION 3.12)

project(pyxtp)

find_package(pybind11 REQUIRED)
pybind11_add_module(pyxtp pyxtp.cc)
    

# DPYXTP_VERSION_INFO is defined in the setup.py file
target_compile_definitions(pyxtp PRIVATE VERSION_INFO=${DPYXTP_VERSION_INFO})