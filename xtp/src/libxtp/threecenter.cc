/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/threecenter.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/symmetric_matrix.h"

namespace votca {
namespace xtp {

void TCMatrix_gwbse::Initialize(Index basissize, Index mmin, Index mmax,
                                Index nmin, Index nmax) {

  // here as storage indices starting from zero
  nmin_ = nmin;
  nmax_ = nmax;
  ntotal_ = nmax - nmin + 1;
  mmin_ = mmin;
  mmax_ = mmax;
  mtotal_ = mmax - mmin + 1;
  auxbasissize_ = basissize;

  // vector has mtotal elements
  // largest object should be allocated in multithread fashion
  matrix_ = std::vector<Eigen::MatrixXd>(mtotal_);
#pragma omp parallel for schedule(dynamic, 4)
  for (Index i = 0; i < mtotal_; i++) {
    matrix_[i] = Eigen::MatrixXd::Zero(ntotal_, auxbasissize_);
  }
}

/*
 * Modify 3-center matrix elements consistent with use of symmetrized
 * Coulomb interaction using either CUDA or Openmp.
 */
void TCMatrix_gwbse::MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix) {
  OpenMP_CUDA gemm;
  gemm.setOperators(matrix_, matrix);
#pragma omp parallel
  {
    Index threadid = OPENMP::getThreadId();
#pragma omp for schedule(dynamic)
    for (Index i = 0; i < msize(); i++) {
      gemm.MultiplyRight(matrix_[i], threadid);
    }
  }
}
/*
 * Fill the 3-center object by looping over shells of GW basis set and
 * calling FillBlock, which calculates all 3-center overlap integrals
 * associated to a particular shell, convoluted with the DFT orbital
 * coefficients
 */
void TCMatrix_gwbse::Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
                          const Eigen::MatrixXd& dft_orbitals) {
  // needed for Rebuild())
  auxbasis_ = &auxbasis;
  dftbasis_ = &dftbasis;
  dft_orbitals_ = &dft_orbitals;

  Fill3cMO(auxbasis, dftbasis, dft_orbitals);

  AOOverlap auxoverlap;
  auxoverlap.Fill(auxbasis);
  AOCoulomb auxcoulomb;
  auxcoulomb.Fill(auxbasis);
  Eigen::MatrixXd inv_sqrt = auxcoulomb.Pseudo_InvSqrt_GWBSE(auxoverlap, 5e-7);
  removedfunctions_ = auxcoulomb.Removedfunctions();
  MultiplyRightWithAuxMatrix(inv_sqrt);

  return;
}

// =============================================================================
// TCMatrix_gwbse::Rotate
//
// Rotates only the n-index (inner rows) of ALL m-slices for the QP window.
//
// In QSGW the self-energy matrix element indices (outer m) stay in the
// DFT-MO basis. Only the construction sum indices (inner n-rows) need to be
// in the QP wavefunction basis. This method applies:
//
//   new_M[m].middleRows(qp_offset_n, qptotal) =
//       U^T * old_M[m].middleRows(qp_offset_n, qptotal)
//
// for ALL m-slices in the full RPA range (because every sigma calculation
// uses every m-slice as a matrix element index and needs updated n-rows).
// Rows outside the QP window remain as DFT-MOs.
// =============================================================================
void TCMatrix_gwbse::Rotate(const Eigen::MatrixXd& U, Index qpmin, Index qpmax) {
  const Index qptotal     = qpmax - qpmin + 1;
  const Index qp_offset_n = qpmin - nmin_;

  assert(qpmin >= nmin_ && qpmax <= nmax_);
  assert(U.rows() == qptotal && U.cols() == qptotal);

  // Rotate the QP-window block of n-rows in every m-slice.
  // new_rows = U^T * old_rows  (U^T left-multiplies the row block)
#pragma omp parallel for schedule(dynamic)
  for (Index m = 0; m < mtotal_; m++) {
    matrix_[m].middleRows(qp_offset_n, qptotal) =
        U.transpose() * matrix_[m].middleRows(qp_offset_n, qptotal);
  }
}

}  // namespace xtp
}  // namespace votca
