/*
 *            Copyright 2009-2026 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_ACCELERATE_LAPACK_H
#define VOTCA_XTP_ACCELERATE_LAPACK_H

namespace votca {
namespace xtp {

/**
 * \brief LAPACK routines of Apple Accelerate, called directly.
 *
 * With Accelerate, Eigen's LAPACK calls go through LAPACKE, and the LAPACKE
 * the build finds (e.g. Homebrew's) sits on a reference LAPACK and BLAS:
 * single threaded and slow. These wrappers call Accelerate's own LAPACK
 * (the interface of macOS 13.3 and later) instead. They live in a file that
 * does not see lapacke.h, so the two LAPACK declarations never meet.
 *
 * Column-major double matrices, leading dimension lda. Each function returns
 * the LAPACK info (0 on success), or kAccelerateUnavailable when this build
 * or this macOS has no Accelerate LAPACK; callers then use their other path.
 */
namespace accelerate {

constexpr int kAccelerateUnavailable = -1000;

/// Whether the routines below can run here (Accelerate build, macOS 13.3+).
bool LapackAvailable();

/// Eigenvalues (ascending, into w) and, for jobz 'V', eigenvectors (into a)
/// of a symmetric matrix, divide and conquer.
int dsyevd(char jobz, char uplo, long n, double* a, long lda, double* w);

/// Cholesky factorisation of a symmetric positive definite matrix.
int dpotrf(char uplo, long n, double* a, long lda);

/// Inverse from the Cholesky factor of dpotrf (one triangle).
int dpotri(char uplo, long n, double* a, long lda);

}  // namespace accelerate
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ACCELERATE_LAPACK_H
