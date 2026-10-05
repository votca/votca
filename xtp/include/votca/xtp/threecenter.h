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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_THREECENTER_H
#define VOTCA_XTP_THREECENTER_H

// Local VOTCA includes
#include "aobasis.h"
#include "eigen.h"
#include "symmetric_matrix.h"

/**
 * \brief Calculates three electron repulsion integrals for GW and DFT.
 *
 *
 *
 */

namespace votca {
namespace xtp {

// due to different requirements for the data format for DFT and GW we have two
// different classes TCMatrix_gwbse and TCMatrix_dft which inherit from TCMatrix
class TCMatrix {

 public:
  TCMatrix() = default;
  virtual ~TCMatrix() = default;
  // the virtual destructor would otherwise suppress moving: moving the DFT
  // tensor (tens of GB) must not fall back to a copy
  TCMatrix(const TCMatrix&) = default;
  TCMatrix(TCMatrix&&) = default;
  TCMatrix& operator=(const TCMatrix&) = default;
  TCMatrix& operator=(TCMatrix&&) = default;
  Index Removedfunctions() const { return removedfunctions_; }
  enum class SpinChannel { Alpha, Beta };

 protected:
  Index removedfunctions_ = 0;
  Eigen::MatrixXd inv_sqrt_;
};

/**
 * \brief RI three-centre tensor of the DFT code,
 * B_P,mu nu = sum_Q V^-1/2_PQ (Q|mu nu).
 *
 * Only the significant basis-function pairs are stored. A shell pair (a,b) is
 * kept if the Schwarz bound of its three-centre integrals,
 * max_P sqrt((P|P)) sqrt((ab|ab)), reaches the pair threshold; the integrals
 * of the dropped pairs are smaller than the threshold. As V^-1/2 only mixes
 * the aux index, the dropped pairs are the same for every aux function.
 * Threshold 0 keeps every pair.
 *
 * Storage: one column of data_ per aux function, holding the kept (nu, mu)
 * with nu <= mu, column mu after column mu-1 (upper triangle, column-major).
 */
class TCMatrix_dft final : public TCMatrix {
 public:
  /// build_columns: pair columns one thread builds at a time (see
  /// kBuildColumns); only changes the memory and work split, not the result.
  void Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
            double pair_threshold = 0.0, Index build_columns = kBuildColumns);

  /// Number of aux functions.
  Index size() const { return Index(data_.cols()); }
  /// Number of basis functions.
  Index basissize() const { return basissize_; }

  /// B_P as a full symmetric matrix, zero for the dropped pairs.
  Eigen::MatrixXd FullMatrix(Index P) const;

  /// Writes B_P into full (basissize x basissize). The entries of dropped
  /// pairs are not written: full must hold zeros there, as it does when the
  /// same buffer is reused for all aux functions.
  void FillFullMatrix(Index P, Eigen::MatrixXd& full) const;

  /// Packs the lower triangle of a matrix into the stored pair layout,
  /// off-diagonal entries counted twice, so that for symmetric D
  /// PackWeighted(D).dot(data column P) = sum_{mu nu} B_P,mu nu D_mu nu.
  /// Like Symmetric_Matrix(full), the upper triangle is not read.
  Eigen::VectorXd PackWeighted(const Eigen::MatrixXd& full) const;

  /// Inverse layout map without weights: symmetric matrix from a packed
  /// vector, zero for the dropped pairs.
  Eigen::MatrixXd Unpack(const Eigen::VectorXd& packed) const;

  /// The stored tensor, one column per aux function.
  const Eigen::MatrixXd& Data() const { return data_; }

  /// Stored pairs (nu <= mu) and all pairs N(N+1)/2.
  Index StoredPairs() const { return Index(data_.rows()); }
  Index AllPairs() const { return basissize_ * (basissize_ + 1) / 2; }

  double MetricSeconds() const { return metric_seconds_; }

 private:
  /// Contiguous piece of column mu: rows [row, row+length) of B, stored from
  /// position offset of a data column on.
  struct Run {
    Index col;
    Index row;
    Index length;
    Index offset;
  };

  /// Kept shells of the pairs (a, b) with b <= a, for each shell a, in
  /// increasing order; a itself is always kept.
  std::vector<std::vector<Index>> SignificantPairs(
      const AOBasis& dftbasis, const Eigen::MatrixXd& metric,
      double pair_threshold) const;
  void SetupLayout(const AOBasis& dftbasis,
                   const std::vector<std::vector<Index>>& kept);

  Eigen::MatrixXd data_;
  std::vector<Run> runs_;
  // column mu starts at col_offset_[mu] in a data column
  std::vector<Index> col_offset_;
  Index basissize_ = 0;
  double metric_seconds_ = 0.0;

  /// Upper bound for the pair columns that one thread builds at a time:
  /// the transient memory of Fill is two such blocks of aux size per thread.
  static constexpr Index kBuildColumns = 512;
};

class TCMatrix_gwbse final : public TCMatrix {
 public:
  // Eigenvalue tolerance of the metric (Pseudo_InvSqrt_GWBSE) folded into
  // the stored integrals. Shared with EnvironmentScreening::Metric, which
  // must build the same T.
  static constexpr double metric_tolerance = 5e-7;

  // returns one level as a constant reference
  const Eigen::MatrixXd& operator[](Index i) const { return matrix_[i]; }

  // returns one level as a reference
  Eigen::MatrixXd& operator[](Index i) { return matrix_[i]; }
  // returns auxbasissize
  Index auxsize() const { return auxbasissize_; }

  Index get_mmin() const { return mmin_; }

  Index get_mmax() const { return mmax_; }

  Index get_nmin() const { return nmin_; }

  Index get_nmax() const { return nmax_; }

  Index msize() const { return mtotal_; }

  Index nsize() const { return ntotal_; }

  void Initialize(Index basissize, Index mmin, Index mmax, Index nmin,
                  Index nmax);

  void Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
            const Eigen::MatrixXd& dft_orbitals);
  // Rebuilds ThreeCenterIntegrals, only works if the original basisobjects
  // still exist
  void Rebuild() { Fill(*auxbasis_, *dftbasis_, *dft_orbitals_); }

  void MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix);

  // The metric that was folded into the stored integrals at Fill time:
  // T = Pseudo_InvSqrt_GWBSE, with T T^T = V^-1 (pseudo-inverse) and
  // T^T V T the projector onto the retained auxiliary functions. After
  // Fill, operator[] holds M = (mn|Q) T, so the bare Coulomb operator
  // exists only implicitly, as M M^T.
  //
  // Kept because anything that must act on the Coulomb interaction in
  // the same space as M needs it. The environment reaction field is the
  // case in point: a kernel B between auxiliary functions, taken as charge
  // densities, enters as M (T^T B T) M^T -- see EnvironmentScreening.
  //
  // IN THE FILL-TIME AUXILIARY BASIS. MultiplyRightWithAuxMatrix, which the
  // PPM and the BSE both call with an eigenvector matrix U, rotates the
  // auxiliary index of M. This matrix is not rotated with it: a kernel
  // built from it has to be brought into the current frame with
  // ToCurrentAuxFrame before it meets operator[]. Refilled, and so reset,
  // by every Fill and Rebuild.
  const Eigen::MatrixXd& InvSqrt() const { return inv_sqrt_; }

  // The auxiliary frame operator[] is currently in. Every
  // MultiplyRightWithAuxMatrix since the last Fill or Rebuild is recorded
  // here, so that operator[] holds M_fill * AuxFrame() with M_fill the
  // integrals as Fill left them. Empty means the identity: nothing has
  // rotated them yet.
  //
  // The bare Coulomb interaction does not care -- it is the identity in
  // every orthonormal frame, which is why nothing tracked this before. A
  // kernel K between the columns of M_fill does care: M_fill K M_fill^T is
  // M K' M^T only with K' = ToCurrentAuxFrame(K).
  const Eigen::MatrixXd& AuxFrame() const { return aux_frame_; }

  // For code that snapshots slices through operator[] and later writes
  // them back, which bypasses the bookkeeping above: restore the frame
  // that was current when the snapshot was taken, alongside the slices.
  void RestoreAuxFrame(const Eigen::MatrixXd& frame) {
    aux_frame_ = frame;
    aux_frame_known_ = true;
  }

  // K, a kernel between the columns of M_fill, as it must be applied to
  // operator[] now: F^-1 K F^-T with F = AuxFrame(), which is F^T K F for
  // the orthogonal eigenvector matrices the PPM and the BSE rotate with.
  // Throws if the frame is not known -- after a non-square
  // MultiplyRightWithAuxMatrix, which is not a change of frame.
  Eigen::MatrixXd ToCurrentAuxFrame(const Eigen::MatrixXd& K) const;

  // Whether AuxFrame() is orthogonal, so that the bare interaction M M^T is
  // what it was at Fill time. False once the index has been dressed.
  bool AuxFrameIsOrthogonal() const;

  /**
   * \brief Dress the auxiliary index with S: M_fill -> M_fill S.
   *
   * For a static environment the QM electrons interact through
   * u = v + v_reac, which is 1 + R in the metric of M, and the screened
   * interaction is W = [u^-1 - chi0]^-1. With S = (1 + R)^(1/2), that is
   * exactly the ordinary RPA and correlation self-energy built from M S
   * in place of M: every place that pairs two M's through the bare
   * interaction pairs them through u instead.
   *
   * S is in the fill-time frame, and is applied in whatever frame the
   * integrals are in now; AuxFrame() records it, so ToCurrentAuxFrame
   * still maps fill-time kernels correctly. Exchange-type sums that must
   * keep the bare v take ToCurrentAuxFrame(identity) as their kernel once
   * the frame is no longer orthogonal (see AuxFrameIsOrthogonal).
   *
   * Throws if already dressed. Fill and Rebuild undress.
   */
  void DressAuxIndex(const Eigen::MatrixXd& S);
  // Back to M_fill (times whatever rotations happened since), bare v.
  void UndressAuxIndex();
  bool Dressed() const { return dressing_.size() > 0; }

  /**
   * \brief Rotate the n-index (construction rows) of Mmn for QSGW.
   *
   * In QSGW the self-energy matrix elements are in the DFT-MO basis
   * (outer m-index unchanged). The internal construction sum over particle/hole
   * states (inner n-rows) must use QP wavefunctions. This method rotates only
   * the n-rows of ALL m-slices within the QP window block:
   *
   *   new_M[m].rows(qp_offset_n : qp_offset_n+qptotal) =
   *       U^T * old_M[m].rows(qp_offset_n : qp_offset_n+qptotal)
   *
   * Rows outside the QP window remain as DFT-MOs (consistent with evGW).
   * The outer m-index and auxiliary index are unchanged.
   *
   * @param U      Rotation matrix (qptotal x qptotal)
   * @param qpmin  First QP level (absolute MO index)
   * @param qpmax  Last  QP level (absolute MO index)
   */
  void Rotate(const Eigen::MatrixXd& U, Index qpmin, Index qpmax);

 private:
  // store vector of matrices
  std::vector<Eigen::MatrixXd> matrix_;

  // band summation indices
  Index mmin_;
  Index mmax_;
  Index nmin_;
  Index nmax_;
  Index ntotal_;
  Index mtotal_;
  Index auxbasissize_;

  Eigen::MatrixXd aux_frame_;  // empty: identity
  bool aux_frame_known_ = true;
  Eigen::MatrixXd dressing_;  // S of DressAuxIndex; empty: bare

  // M_fill -> M_fill A, whatever frame the integrals are in.
  void MultiplyInFillFrame(const Eigen::MatrixXd& A);

  const AOBasis* auxbasis_ = nullptr;
  const AOBasis* dftbasis_ = nullptr;
  const Eigen::MatrixXd* dft_orbitals_ = nullptr;

  void Fill3cMO(const AOBasis& auxbasis, const AOBasis& dftbasis,
                const Eigen::MatrixXd& dft_orbitals);
};

struct TCMatrix_gwbse_spin {
  TCMatrix_gwbse alpha;
  TCMatrix_gwbse beta;

  TCMatrix_gwbse& operator[](TCMatrix::SpinChannel spin) {
    return (spin == TCMatrix::SpinChannel::Alpha) ? alpha : beta;
  }

  const TCMatrix_gwbse& operator[](TCMatrix::SpinChannel spin) const {
    return (spin == TCMatrix::SpinChannel::Alpha) ? alpha : beta;
  }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_THREECENTER_H
