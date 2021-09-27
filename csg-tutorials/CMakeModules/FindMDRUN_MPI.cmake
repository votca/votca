# - this module looks for mpi-enable mdrun
#
# Once done this will define
#
#  MDRUN_MPI_FOUND      - system has mdrun_mpi
#  MDRUN_MPI_EXECUTABLE - the mprun_mpi executable
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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

INCLUDE(FindCygwin)

FIND_PROGRAM(MDRUN_MPI_EXECUTABLE
  NAMES 
  mdrun_mpi_d
  mdrun_openmpi_d #fedora
  mdrun_mpi.openmpi_d #debian
  mdrun_mpi
  mdrun_openmpi #fedora
  mdrun_mpi.openmpi #debian
  PATHS
  ${CYGWIN_INSTALL_PATH}/bin
)

# handle the QUIETLY and REQUIRED arguments and set MDRUN_MPI_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MDRUN_MPI DEFAULT_MSG MDRUN_MPI_EXECUTABLE)

MARK_AS_ADVANCED( MDRUN_MPI_EXECUTABLE )

