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
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>

#include "votca/xtp/aomatrix.h"
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
 * Coulomb interaction.
 */
void TCMatrix_gwbse::MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix) {
#pragma omp parallel for schedule(dynamic)
  for (Index i = 0; i < msize(); i++) {
    matrix_[i] *= matrix;
  }
  if (matrix.rows() != matrix.cols()) {
    aux_frame_known_ = false;  // a projection, not a change of frame
    aux_frame_.resize(0, 0);
  } else if (aux_frame_known_) {
    aux_frame_ = (aux_frame_.size() == 0) ? matrix : aux_frame_ * matrix;
  }
}

void TCMatrix_gwbse::MultiplyRightWithAuxMatrixLeading(const Eigen::MatrixXd& U,
                                                       Index full_slices,
                                                       Index lead) {
  if (U.rows() != U.cols()) {
    MultiplyRightWithAuxMatrix(U);
    return;
  }
  const Index lead_rows = std::clamp<Index>(lead, 0, nsize());
  const Index lead_slices = std::clamp<Index>(lead, 0, msize());
  full_slices = std::clamp<Index>(full_slices, 0, lead_slices);
#pragma omp parallel for schedule(dynamic)
  for (Index i = 0; i < lead_slices; i++) {
    const Index rows = (i < full_slices) ? nsize() : lead_rows;
    Eigen::MatrixXd rotated = matrix_[i].topRows(rows) * U;
    matrix_[i].topRows(rows) = rotated;
  }
  if (aux_frame_known_) {
    aux_frame_ = (aux_frame_.size() == 0) ? U : aux_frame_ * U;
  }
  if (full_slices < msize() && (lead_slices < msize() || lead_rows < nsize())) {
    partially_rotated_ = true;
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
  if (partially_rotated_) {
    throw std::runtime_error(
        "TCMatrix_gwbse: cannot dress the auxiliary index after only part of "
        "the tensor was rotated (MultiplyRightWithAuxMatrixLeading). Rebuild "
        "first.");
  }
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
  partially_rotated_ = false;
  dressing_.resize(0, 0);

  return;
}

void TCMatrix_gwbse::RotateOrbitals(const Eigen::MatrixXd& U, Index first,
                                    Index last) {
  const Index window = last - first + 1;
  assert(first >= nmin_ && last <= nmax_);
  assert(first >= mmin_ && last <= mmax_);
  assert(U.rows() == window && U.cols() == window);
  const Index row0 = first - nmin_;
  const Index slice0 = first - mmin_;
  const Index rows = matrix_.empty() ? 0 : Index(matrix_[0].rows());
  const Index aux = matrix_.empty() ? 0 : Index(matrix_[0].cols());

  // second index, every slice
  const Eigen::MatrixXd Ut = U.transpose();
#pragma omp parallel for schedule(dynamic)
  for (Index m = 0; m < Index(matrix_.size()); m++) {
    auto block = matrix_[m].middleRows(row0, window);
    block = Ut * block;
  }

  // first index: the slices of the window as columns of one matrix, in
  // chunks of auxiliary columns (contiguous in each slice) of about 64 MB
  const Index chunk = std::max<Index>(
      1, std::min<Index>(aux, Index(8e6 / double(rows * window + 1))));
  Eigen::MatrixXd X;
  Eigen::MatrixXd Y;
  for (Index c0 = 0; c0 < aux; c0 += chunk) {
    const Index nc = std::min(chunk, aux - c0);
    const Index len = rows * nc;
    X.resize(len, window);
#pragma omp parallel for
    for (Index k = 0; k < window; k++) {
      X.col(k) = Eigen::Map<const Eigen::VectorXd>(
          matrix_[slice0 + k].col(c0).data(), len);
    }
    Y.noalias() = X * U;
#pragma omp parallel for
    for (Index i = 0; i < window; i++) {
      Eigen::Map<Eigen::VectorXd>(matrix_[slice0 + i].col(c0).data(), len) =
          Y.col(i);
    }
  }
}

void TCMatrix_dft::SetupLayout(const AOBasis& dftbasis,
                               const std::vector<std::vector<Index>>& kept) {
  basissize_ = dftbasis.AOBasisSize();
  const std::vector<Index> shell2bf = dftbasis.getMapToBasisFunctions();
  runs_.clear();
  col_offset_.assign(basissize_, 0);
  Index offset = 0;
  for (Index a = 0; a < Index(kept.size()); ++a) {
    const Index start = shell2bf[a];
    const Index size = dftbasis.getShell(a).getNumFunc();
    for (Index i = 0; i < size; ++i) {
      const Index mu = start + i;
      col_offset_[mu] = offset;
      // merge neighbouring kept shells into one run of rows
      Index k = 0;
      while (k < Index(kept[a].size())) {
        Index b_first = kept[a][k];
        Index b_last = b_first;
        while (k + 1 < Index(kept[a].size()) && kept[a][k + 1] == b_last + 1) {
          ++k;
          b_last = kept[a][k];
        }
        ++k;
        const Index row = shell2bf[b_first];
        const Index end = std::min(
            shell2bf[b_last] + dftbasis.getShell(b_last).getNumFunc(), mu + 1);
        runs_.push_back({mu, row, end - row, offset});
        offset += end - row;
      }
    }
  }
}

void TCMatrix_dft::FillFullMatrix(Index P, Eigen::MatrixXd& full) const {
  assert(full.rows() == basissize_ && full.cols() == basissize_);
  const double* src = data_.col(P).data();
  for (const Run& r : runs_) {
    std::copy(src + r.offset, src + r.offset + r.length,
              full.data() + r.col * basissize_ + r.row);
  }
  full.triangularView<Eigen::StrictlyLower>() = full.transpose();
}

Eigen::MatrixXd TCMatrix_dft::FullMatrix(Index P) const {
  Eigen::MatrixXd full = Eigen::MatrixXd::Zero(basissize_, basissize_);
  FillFullMatrix(P, full);
  return full;
}

Eigen::VectorXd TCMatrix_dft::PackWeighted(const Eigen::MatrixXd& full) const {
  assert(full.rows() == basissize_ && full.cols() == basissize_);
  Eigen::VectorXd packed(data_.rows());
  for (const Run& r : runs_) {
    for (Index k = 0; k < r.length; ++k) {
      const Index row = r.row + k;
      packed(r.offset + k) =
          (row == r.col) ? full(r.col, row) : 2.0 * full(r.col, row);
    }
  }
  return packed;
}

Eigen::MatrixXd TCMatrix_dft::Unpack(const Eigen::VectorXd& packed) const {
  assert(packed.size() == data_.rows());
  Eigen::MatrixXd full = Eigen::MatrixXd::Zero(basissize_, basissize_);
  for (const Run& r : runs_) {
    full.block(r.row, r.col, r.length, 1) = packed.segment(r.offset, r.length);
  }
  full.triangularView<Eigen::StrictlyLower>() = full.transpose();
  return full;
}

}  // namespace xtp
}  // namespace votca
