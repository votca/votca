# Find sphinx-build
#
# This will define
#
# SPHINX_FOUND - Sphinx documentation builder is installed
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
find_program(SPHINX_EXECUTABLE NAMES sphinx-build DOC "Sphinx documentation generation tool (http://www.sphinx-doc.org/)")

if(SPHINX_EXECUTABLE)
execute_process(COMMAND ${SPHINX_EXECUTABLE} --version
  OUTPUT_VARIABLE sphinx_version
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(sphinx_version MATCHES "^sphinx-build .*")
    string(REGEX
           REPLACE "sphinx-build ([.0-9]+).*"
                   "\\1"
                   SPHINX_VERSION
                   "${sphinx_version}")
  else()
    set(SPHINX_VERSION 0.0)
  endif()
else()
  set(SPHINX_VERSION 0.0)
endif()
message(status "SPHINX_VERSION:${SPHINX_VERSION}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPHINX REQUIRED_VARS SPHINX_EXECUTABLE VERSION_VAR SPHINX_VERSION)
