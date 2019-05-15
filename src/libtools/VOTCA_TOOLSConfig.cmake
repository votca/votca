include(CMakeFindDependencyMacro)
find_dependency(Eigen3)
find_dependency(Boost 1.53.0 REQUIRED COMPONENTS program_options)
include("${CMAKE_CURRENT_LIST_DIR}/VOTCA_TOOLS_Targets.cmake")
