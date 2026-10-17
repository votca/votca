# Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
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

# Apple's Accelerate framework: BLAS for Eigen (EIGEN_USE_BLAS) and, from
# macOS 13.3 on, LAPACK through accelerate_lapack.h (dsyevd, dpotrf,
# dpotri for the large eigensolvers and Cholesky factorisations).
#
# No LAPACKE: Accelerate has no LAPACKE interface, and a Homebrew lapacke
# links the reference LAPACK, whose BLAS calls are bound (two-level
# namespace) to the reference BLAS it was built against, not to
# Accelerate. Every Eigen decomposition routed there (EIGEN_USE_LAPACKE)
# ran single-threaded on unoptimised BLAS. Eigen's own decompositions use
# the overloaded BLAS-3 kernels instead.
#
# Defines APPLE_ACCELERATE_FOUND and the imported target
# APPLE_ACCELERATE::Accelerate.

find_library(ACCELERATE Accelerate)
mark_as_advanced(ACCELERATE)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(APPLE_ACCELERATE DEFAULT_MSG ACCELERATE)

if(APPLE_ACCELERATE_FOUND AND NOT TARGET APPLE_ACCELERATE::Accelerate)
  add_library(APPLE_ACCELERATE::Accelerate INTERFACE IMPORTED)
  set_target_properties(APPLE_ACCELERATE::Accelerate PROPERTIES
    INTERFACE_LINK_LIBRARIES "${ACCELERATE}")
endif()
