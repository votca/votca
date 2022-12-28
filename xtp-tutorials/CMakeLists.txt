cmake_minimum_required(VERSION 3.13)

project(xtp-tutorials LANGUAGES NONE)

set(PROJECT_VERSION "2023-dev")

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ DESTINATION ${CMAKE_INSTALL_DATADIR}/votca/xtp-tutorials USE_SOURCE_PERMISSIONS
	PATTERN ".git*" EXCLUDE
	PATTERN "CMake*" EXCLUDE PATTERN "cmake_install.cmake" EXCLUDE
	PATTERN "install_manifest.txt" EXCLUDE PATTERN "*~" EXCLUDE)

add_subdirectory(share/doc)
