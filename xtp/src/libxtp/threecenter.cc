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
// Standard includes
#include <stdexcept>
#include <string>

#include "votca/xtp/aomatrix.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/symmetric_matrix.h"
#include "votca/xtp/threecenter.h"

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
  if (matrix.rows() != matrix.cols()) {
    aux_frame_known_ = false;  // a projection, not a change of frame
    aux_frame_.resize(0, 0);
  } else if (aux_frame_known_) {
    aux_frame_ = (aux_frame_.size() == 0) ? matrix : aux_frame_ * matrix;
  }
}

bool TCMatrix_gwbse::AuxFrameIsOrthogonal() const {
  if (!aux_frame_known_) {
    return false;
  }
  if (aux_frame_.size() == 0) {
    return true;
  }
  const Index n = aux_frame_.rows();
  return (aux_frame_.transpose() * aux_frame_ - Eigen::MatrixXd::Identity(n, n))
             .cwiseAbs()
             .maxCoeff() < 1e-10;
}

void TCMatrix_gwbse::MultiplyInFillFrame(const Eigen::MatrixXd& A) {
  if (!aux_frame_known_) {
    throw std::runtime_error(
        "TCMatrix_gwbse: cannot dress the auxiliary index after a non-square "
        "MultiplyRightWithAuxMatrix. Rebuild first.");
  }
  if (aux_frame_.size() == 0) {
    MultiplyRightWithAuxMatrix(A);
    return;
  }
  // M_fill F X = M_fill A F  =>  X = F^-1 A F. The frame becomes A F.
  const Eigen::PartialPivLU<Eigen::MatrixXd> lu(aux_frame_);
  MultiplyRightWithAuxMatrix(lu.solve(A * aux_frame_));
}

void TCMatrix_gwbse::DressAuxIndex(const Eigen::MatrixXd& S) {
  if (Dressed()) {
    throw std::runtime_error(
        "TCMatrix_gwbse::DressAuxIndex: already dressed. Undress first.");
  }
  if (S.rows() != auxbasissize_ || S.cols() != auxbasissize_) {
    throw std::runtime_error(
        "TCMatrix_gwbse::DressAuxIndex: S is " + std::to_string(S.rows()) +
        "x" + std::to_string(S.cols()) + ", auxiliary basis has " +
        std::to_string(auxbasissize_) + " functions.");
  }
  MultiplyInFillFrame(S);
  dressing_ = S;
}

void TCMatrix_gwbse::UndressAuxIndex() {
  if (!Dressed()) {
    return;
  }
  MultiplyInFillFrame(dressing_.inverse());
  dressing_.resize(0, 0);
}

Eigen::MatrixXd TCMatrix_gwbse::ToCurrentAuxFrame(
    const Eigen::MatrixXd& K) const {
  if (!aux_frame_known_) {
    throw std::runtime_error(
        "TCMatrix_gwbse::ToCurrentAuxFrame: the auxiliary index of the "
        "three-centre integrals was multiplied by a non-square matrix, so "
        "there is no frame to bring the kernel into. Apply it before that, "
        "or Rebuild.");
  }
  if (aux_frame_.size() == 0) {
    return K;
  }
  if (K.rows() != aux_frame_.rows() || K.cols() != aux_frame_.rows()) {
    throw std::runtime_error(
        "TCMatrix_gwbse::ToCurrentAuxFrame: kernel is " +
        std::to_string(K.rows()) + "x" + std::to_string(K.cols()) +
        ", auxiliary frame is " + std::to_string(aux_frame_.rows()) + "x" +
        std::to_string(aux_frame_.cols()) + ".");
  }
  const Index n = aux_frame_.rows();
  const double off_orthogonal =
      (aux_frame_.transpose() * aux_frame_ - Eigen::MatrixXd::Identity(n, n))
          .cwiseAbs()
          .maxCoeff();
  Eigen::MatrixXd result;
  if (off_orthogonal < 1e-10) {
    result = aux_frame_.transpose() * K * aux_frame_;
  } else {
    const Eigen::PartialPivLU<Eigen::MatrixXd> lu(aux_frame_);
    const Eigen::MatrixXd X = lu.solve(K);  // F^-1 K
    result = lu.solve(Eigen::MatrixXd(X.transpose())).transpose();
  }
  return result;
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
  inv_sqrt_ = auxcoulomb.Pseudo_InvSqrt_GWBSE(auxoverlap, metric_tolerance);
  removedfunctions_ = auxcoulomb.Removedfunctions();
  MultiplyRightWithAuxMatrix(inv_sqrt_);
  // What Fill leaves behind is the reference frame by definition.
  aux_frame_.resize(0, 0);
  aux_frame_known_ = true;
  dressing_.resize(0, 0);

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
void TCMatrix_gwbse::Rotate(const Eigen::MatrixXd& U, Index qpmin,
                            Index qpmax) {
  const Index qptotal = qpmax - qpmin + 1;
  const Index qp_offset_n = qpmin - nmin_;  // row offset in n-storage
  const Index qp_offset_m = qpmin - mmin_;  // slice offset in m-storage

  assert(qpmin >= nmin_ && qpmax <= nmax_);
  assert(qpmin >= mmin_ && qpmax <= mmax_);
  assert(U.rows() == qptotal && U.cols() == qptotal);

  // Rotate the n-rows of ONLY the QP-window m-slices [qpmin, qpmax].
  // Slices outside this range (e.g. core levels below qpmin, or high virtuals
  // above qpmax) are always DFT-MOs and must NOT be touched -- they are not
  // saved/restored by Mmn_orig in gw.cc and would accumulate drift if rotated.
  // Only the QP-window n-row block [qp_offset_n, qp_offset_n+qptotal) is
  // rotated; rows outside remain as DFT-MOs (consistent with evGW treatment).
  // new_rows = U^T * old_rows
#pragma omp parallel for schedule(dynamic)
  for (Index m = 0; m < qptotal; m++) {
    matrix_[m + qp_offset_m].middleRows(qp_offset_n, qptotal) =
        U.transpose() *
        matrix_[m + qp_offset_m].middleRows(qp_offset_n, qptotal);
  }
}

}  // namespace xtp
}  // namespace votca
