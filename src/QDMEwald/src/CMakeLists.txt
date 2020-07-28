add_library(QDMEwald INTERFACE)
target_include_directories(QDMEwald INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(QDMEwald INTERFACE Kokkos::kokkos)

add_executable(dummy dummy.cpp)
target_link_libraries(dummy QDMEwald)
