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
#include "logger.h"
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
  Index Removedfunctions() const { return _removedfunctions; }

 protected:
  Index _removedfunctions = 0;
  Eigen::MatrixXd _inv_sqrt;
};

class TCMatrix_dft final : public TCMatrix {
 public:
  void Fill(const AOBasis& auxbasis, const AOBasis& dftbasis);

  Index size() const { return Index(_matrix.size()); }

  Symmetric_Matrix& operator[](Index i) { return _matrix[i]; }

  const Symmetric_Matrix& operator[](Index i) const { return _matrix[i]; }

 private:
  std::vector<Symmetric_Matrix> _matrix;

  void FillBlock(std::vector<Eigen::MatrixXd>& block, Index shellindex,
                 const AOBasis& dftbasis, const AOBasis& auxbasis);
};

class TCMatrix_gwbse final : public TCMatrix {
 public:
  // returns one level as a constant reference
  const Eigen::MatrixXd& operator[](Index i) const { return _matrix[i]; }

  // returns one level as a reference
  Eigen::MatrixXd& operator[](Index i) { return _matrix[i]; }
  // returns auxbasissize
  Index auxsize() const { return _auxbasissize; }

  Index get_mmin() const { return _mmin; }

  Index get_mmax() const { return _mmax; }

  Index get_nmin() const { return _nmin; }

  Index get_nmax() const { return _nmax; }

  Index msize() const { return _mtotal; }

  Index nsize() const { return _ntotal; }

  void Initialize(Index basissize, Index mmin, Index mmax, Index nmin,
                  Index nmax);

  void Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
            const Eigen::MatrixXd& dft_orbitals);
  // Rebuilds ThreeCenterIntegrals, only works if the original basisobjects
  // still exist
  void Rebuild() { Fill(*_auxbasis, *_dftbasis, *_dft_orbitals); }

  void MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix);

 private:
  // store vector of matrices
  std::vector<Eigen::MatrixXd> _matrix;

  // band summation indices
  Index _mmin;
  Index _mmax;
  Index _nmin;
  Index _nmax;
  Index _ntotal;
  Index _mtotal;
  Index _auxbasissize;

  const AOBasis* _auxbasis = nullptr;
  const AOBasis* _dftbasis = nullptr;
  const Eigen::MatrixXd* _dft_orbitals = nullptr;

  void Fill3cMO(const AOBasis& auxbasis, const AOBasis& dftbasis,
                const Eigen::MatrixXd& dft_orbitals);

  std::vector<Eigen::MatrixXd> ComputeAO3cBlock(const libint2::Shell& auxshell,
                                                const AOBasis& dftbasis,
                                                libint2::Engine& engine) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_THREECENTER_H
