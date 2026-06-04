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
#include "classicalsegment.h"
#include "eigen.h"
#include "qmmolecule.h"
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
  virtual ~TCMatrix() = default;
  Index Removedfunctions() const { return removedfunctions_; }
  enum class SpinChannel { Alpha, Beta };

 protected:
  Index removedfunctions_ = 0;
  Eigen::MatrixXd inv_sqrt_;
};

class TCMatrix_dft final : public TCMatrix {
 public:
  void Fill(const AOBasis& auxbasis, const AOBasis& dftbasis);

  Index size() const { return Index(matrix_.size()); }

  Symmetric_Matrix& operator[](Index i) { return matrix_[i]; }

  const Symmetric_Matrix& operator[](Index i) const { return matrix_[i]; }

 private:
  std::vector<Symmetric_Matrix> matrix_;

  void FillBlock(std::vector<Eigen::MatrixXd>& block, Index shellindex,
                 const AOBasis& dftbasis, const AOBasis& auxbasis);
};

class TCMatrix_gwbse final : public TCMatrix {
 public:
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

  // Overload that additionally computes and adds the MM reaction-field
  // correction to the 2-centre Coulomb matrix before forming V^{-1/2}.
  // polar_segs must contain Thole-parametrised PolarSites; their
  // induced-dipole fields are used as scratch and are reset on return.
  // grid_quality is passed to Vxc_Grid::GridSetup; "xfine" is recommended.
  void Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
            const Eigen::MatrixXd& dft_orbitals, const QMMolecule& qmatoms,
            std::vector<PolarSegment>& polar_segs,
            const std::string& grid_quality = "xfine");

  // Rebuilds ThreeCenterIntegrals, only works if the original basisobjects
  // still exist
  void Rebuild() { Fill(*auxbasis_, *dftbasis_, *dft_orbitals_); }

  void MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix);

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
