# Find jupyter nbconvert
#
# This will define
#
# JUPYTER_NBCONVERT_FOUND - Jupyter and nbsphinx python modules are installed
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

find_program(JUPYTER_EXECUTABLE NAMES jupyter DOC "Interactive computing environment (https://jupyter.org/)")
find_package_handle_standard_args(JUPYTER REQUIRED_VARS JUPYTER_EXECUTABLE)

set(JUPYTER_NBCONVERT_FOUND FALSE)
if(JUPYTER_FOUND)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import cma"
    RESULT_VARIABLE AVAILABLE_NBCONVERT)
  if(AVAILABLE_NBCONVERT EQUAL 0)
    set(JUPYTER_NBCONVERT_FOUND TRUE)
  endif()
endif()

find_package_handle_standard_args(JUPYTER_NBCONVERT DEFAULT_MSG JUPYTER_NBCONVERT_FOUND)
