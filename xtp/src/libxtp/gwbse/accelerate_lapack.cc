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

// Deliberately no Eigen or VOTCA eigen.h here: those pull in lapacke.h,
// whose LAPACK prototypes clash with Accelerate's.
#include <votca/tools/votca_tools_config.h>

#include "votca/xtp/accelerate_lapack.h"

// The new LAPACK interface of Accelerate needs the macOS 13.3 SDK; which
// macOS the program then runs on is checked at run time.
#if defined(__APPLE__) && defined(__clang__) && defined(APPLE_ACCELERATE_FOUND)
#include <Availability.h>
#if defined(__MAC_OS_X_VERSION_MAX_ALLOWED) && \
    __MAC_OS_X_VERSION_MAX_ALLOWED >= 130300
#define VOTCA_XTP_ACCELERATE_LAPACK
#endif
#endif

#ifdef VOTCA_XTP_ACCELERATE_LAPACK
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#include <algorithm>
#include <vector>
#endif

namespace votca {
namespace xtp {
namespace accelerate {

#ifdef VOTCA_XTP_ACCELERATE_LAPACK

bool LapackAvailable() {
  if (__builtin_available(macOS 13.3, *)) {
    return true;
  }
  return false;
}

int dsyevd(char jobz, char uplo, long n, double* a, long lda, double* w) {
  if (__builtin_available(macOS 13.3, *)) {
    const __LAPACK_int N = static_cast<__LAPACK_int>(n);
    const __LAPACK_int LDA = static_cast<__LAPACK_int>(std::max(1L, lda));
    __LAPACK_int info = 0;
    // workspace query
    double work_query = 0.0;
    __LAPACK_int iwork_query = 0;
    __LAPACK_int lwork = -1;
    __LAPACK_int liwork = -1;
    dsyevd_(&jobz, &uplo, &N, a, &LDA, w, &work_query, &lwork, &iwork_query,
            &liwork, &info);
    if (info != 0) {
      return static_cast<int>(info);
    }
    lwork = static_cast<__LAPACK_int>(work_query);
    liwork = iwork_query;
    std::vector<double> work(
        static_cast<std::size_t>(std::max<__LAPACK_int>(1, lwork)));
    std::vector<__LAPACK_int> iwork(
        static_cast<std::size_t>(std::max<__LAPACK_int>(1, liwork)));
    dsyevd_(&jobz, &uplo, &N, a, &LDA, w, work.data(), &lwork, iwork.data(),
            &liwork, &info);
    return static_cast<int>(info);
  }
  return kAccelerateUnavailable;
}

int dpotrf(char uplo, long n, double* a, long lda) {
  if (__builtin_available(macOS 13.3, *)) {
    const __LAPACK_int N = static_cast<__LAPACK_int>(n);
    const __LAPACK_int LDA = static_cast<__LAPACK_int>(std::max(1L, lda));
    __LAPACK_int info = 0;
    dpotrf_(&uplo, &N, a, &LDA, &info);
    return static_cast<int>(info);
  }
  return kAccelerateUnavailable;
}

int dpotri(char uplo, long n, double* a, long lda) {
  if (__builtin_available(macOS 13.3, *)) {
    const __LAPACK_int N = static_cast<__LAPACK_int>(n);
    const __LAPACK_int LDA = static_cast<__LAPACK_int>(std::max(1L, lda));
    __LAPACK_int info = 0;
    dpotri_(&uplo, &N, a, &LDA, &info);
    return static_cast<int>(info);
  }
  return kAccelerateUnavailable;
}

#else

bool LapackAvailable() { return false; }
int dsyevd(char, char, long, double*, long, double*) {
  return kAccelerateUnavailable;
}
int dpotrf(char, long, double*, long) { return kAccelerateUnavailable; }
int dpotri(char, long, double*, long) { return kAccelerateUnavailable; }

#endif

}  // namespace accelerate
}  // namespace xtp
}  // namespace votca
