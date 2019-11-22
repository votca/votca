# - this module looks for ESPRESSO
#
# Once done this will define
#
#  ESPRESSO_FOUND      - system has lmp
#
# Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

find_package(PythonInterp 3)
set(IMPORT_ESPRESSO_SUCCESS FALSE)
if(PythonInterp_FOUND)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import espressomd"
    RESULT_VARIABLE IMPORT_ESPRESSO)
  if(IMPORT_ESPRESSO EQUAL 0)
    set(IMPORT_ESPRESSO_SUCCESS TRUE)
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set LMP_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ESPRESSO DEFAULT_MSG IMPORT_ESPRESSO_SUCCESS)
