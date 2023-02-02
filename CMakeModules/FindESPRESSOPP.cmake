# - this module looks for ESPRESSOPP
#
# Once done this will define
#
#  ESPRESSO_FOUND      - system has espressomd
#
# Copyright 2009-2022 The VOTCA Development Team (http://www.votca.org)
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

set(IMPORT_ESPRESSOPP_SUCCESS FALSE)
if(Python_EXECUTABLE)
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import espressopp"
    RESULT_VARIABLE IMPORT_ESPRESSOPP)
  if(IMPORT_ESPRESSOPP EQUAL 0)
    set(IMPORT_ESPRESSOPP_SUCCESS TRUE)
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set LMP_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ESPRESSOPP DEFAULT_MSG IMPORT_ESPRESSOPP_SUCCESS)
