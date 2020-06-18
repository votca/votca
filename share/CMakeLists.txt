add_subdirectory(xtp/xml)
add_subdirectory(man)
add_subdirectory(basis_sets)
add_subdirectory(ecps)
add_subdirectory(benchmark)
add_subdirectory(doc)

file(COPY xtp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY basis_sets DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/xtp)
