# - this module looks for ESPRESSO
#
# Once done this will define
#
#  ESPRESSO_FOUND      - system has espressomd
#
# Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if(CMAKE_VERSION VERSION_LESS 3.12)
  find_package(PythonInterp 3)
  if(PYTHONINTERP_FOUND)
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
  endif()
else()
  find_package(Python COMPONENTS Interpreter)
endif()

set(IMPORT_ESPRESSO_SUCCESS FALSE)
if(Python_EXECUTABLE)
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import espressomd"
    RESULT_VARIABLE IMPORT_ESPRESSO)
  if(IMPORT_ESPRESSO EQUAL 0)
    set(IMPORT_ESPRESSO_SUCCESS TRUE)
  endif()
endif()

if(IMPORT_ESPRESSO_SUCCESS)
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import espressomd; exit(not espressomd.has_features('H5MD'))"
    RESULT_VARIABLE ESPRESSO_HAS_H5MD)
  if(ESPRESSO_HAS_H5MD EQUAL 0)
    set(ESPRESSO_HAS_H5MD_SUCCESS TRUE)
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set LMP_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ESPRESSO DEFAULT_MSG IMPORT_ESPRESSO_SUCCESS ESPRESSO_HAS_H5MD_SUCCESS)
