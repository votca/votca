/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#ifndef __XTP_THREECENTER__H
#define __XTP_THREECENTER__H

#include <array>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/symmetric_matrix.h>

#ifdef USE_CUDA
#include "cudapipeline.h"
#endif

/**
 * \brief Calculates three electron overlap integrals for GW and DFT.
 *
 *
 *
 */

namespace votca {
namespace xtp {

class AOShell;
class AOBasis;

// due to different requirements for the data format for DFT and GW we have two
// different classes TCMatrix_gwbse and TCMatrix_dft which inherit from TCMatrix
class TCMatrix {

 public:
  Index Removedfunctions() const { return _removedfunctions; }

 protected:
  Index _removedfunctions = 0;
  Eigen::MatrixXd _inv_sqrt;

  bool FillThreeCenterRepBlock(Eigen::Tensor<double, 3>& threec_block,
                               const AOShell& shell, const AOShell& shell_row,
                               const AOShell& shell_col) const;
};

class TCMatrix_dft : public TCMatrix {
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

class TCMatrix_gwbse : public TCMatrix {
 public:
  TCMatrix_gwbse(Logger& log) : _log{log} {};

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

  void MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& AuxMatrix);

 private:
  // store vector of matrices
  std::vector<Eigen::MatrixXd> _matrix;

  // Logger
  Logger& _log;

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

  std::vector<Eigen::MatrixXd> FillBlock(
      const std::vector<Eigen::MatrixXd>& symmstorage,
      const Eigen::MatrixXd& dft_orbitals) const;

  void MultiplyRightWithAuxMatrixOpenMP(const Eigen::MatrixXd& AuxMatrix);

  void FillAllBlocksOpenMP(const AOBasis& gwbasis, const AOBasis& dftbasis,
                           const Eigen::MatrixXd& dft_orbitals);

  std::vector<Eigen::MatrixXd> ComputeSymmStorage(
      const AOShell& auxshell, const AOBasis& dftbasis,
      const Eigen::MatrixXd& dft_orbitals) const;

#if defined(USE_CUDA)
  std::array<CudaMatrix, 2> SendDFTMatricesToGPU(
      const Eigen::MatrixXd& dft_orbitals, const CudaPipeline& cuda_pip) const;

  std::array<CudaMatrix, 3> CreateIntermediateCudaMatrices(
      Index basissize, const CudaPipeline& cuda_pip) const;

  void FillAllBlocksCuda(const AOBasis& gwbasis, const AOBasis& dftbasis,
                         const Eigen::MatrixXd& dft_orbitals);

  void MultiplyRightWithAuxMatrixCuda(const Eigen::MatrixXd& matrix);

  std::vector<Eigen::MatrixXd> FillBlockCUDA(
      const std::vector<Eigen::MatrixXd>& symmstorage,
      const std::array<CudaMatrix, 2>& cuda_matrices,
      std::array<CudaMatrix, 3>& cuda_inter_matrices,
      const CudaPipeline& cuda_pip) const;

#endif
};

}  // namespace xtp
}  // namespace votca

#endif /* AOMATRIX_H */
