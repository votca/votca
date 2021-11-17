# - Find libint2
# Find libint2 through pkg-config
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

find_package(PkgConfig)
if(COMMAND set_package_properties)
  set_package_properties(PkgConfig PROPERTIES TYPE RECOMMENDED PURPOSE "Used to detect libint2 package")
endif()

pkg_check_modules(PC_LIBINT2 REQUIRED IMPORTED_TARGET libint2>=2.6)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(libint2 REQUIRED_VARS PC_LIBINT2_INCLUDE_DIRS PC_LIBINT2_LIBRARIES VERSION_VAR PC_LIBINT2_VERSION) 

if(TARGET PkgConfig::PC_LIBINT2 AND NOT Libint2::int2) 
  set_target_properties(PkgConfig::PC_LIBINT2 PROPERTIES IMPORTED_GLOBAL TRUE)
  add_library(Libint2::int2 ALIAS PkgConfig::PC_LIBINT2)  
endif()
