# - this module looks for Rdkit
#
# Once done this will define
#
#  Rdkit_FOUND      - system has lmp
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

find_package(Python 3 COMPONENTS Interpreter)
set(IMPORT_Rdkit_SUCCESS FALSE)
if(Python_EXECUTABLE)
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import rdkit"
    RESULT_VARIABLE IMPORT_Rdkit)
  if(IMPORT_Rdkit EQUAL 0)
    set(IMPORT_Rdkit_SUCCESS TRUE)
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set Rdkit_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Rdkit DEFAULT_MSG IMPORT_Rdkit_SUCCESS)
