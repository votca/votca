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

#ifndef __XTP_THREECENTER__H
#define __XTP_THREECENTER__H

#include <cstddef>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/multiarray.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/symmetric_matrix.h>

/**
 * \brief Calculates three electron overlap integrals for GW and DFT.
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
  int Removedfunctions() const { return _removedfunctions; }

 protected:
  int _removedfunctions = 0;
  Eigen::MatrixXd _inv_sqrt;

  bool FillThreeCenterRepBlock(tensor3d& threec_block, const AOShell& shell,
                               const AOShell& shell_row,
                               const AOShell& shell_col);
};

class TCMatrix_dft : public TCMatrix {
 public:
  void Fill(const AOBasis& auxbasis, const AOBasis& dftbasis);

  int size() const { return _matrix.size(); }

  Symmetric_Matrix& operator[](int i) { return _matrix[i]; }

  const Symmetric_Matrix& operator[](int i) const { return _matrix[i]; }

 private:
  std::vector<Symmetric_Matrix> _matrix;

  void FillBlock(std::vector<Eigen::MatrixXd>& block, int shellindex,
                 const AOBasis& dftbasis, const AOBasis& auxbasis);
};

class TCMatrix_gwbse : public TCMatrix {
 public:
  // returns one level as a constant reference
  const Eigen::MatrixXd& operator[](int i) const { return _matrix[i]; }

  // returns one level as a reference
  Eigen::MatrixXd& operator[](int i) { return _matrix[i]; }
  // returns auxbasissize
  int auxsize() const { return _basissize; }

  int get_mmin() const { return _mmin; }

  int get_mmax() const { return _mmax; }

  int get_nmin() const { return _nmin; }

  int get_nmax() const { return _nmax; }

  int msize() const { return _mtotal; }

  int nsize() const { return _ntotal; }

  void Initialize(int basissize, int mmin, int mmax, int nmin, int nmax);

  void Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
            const Eigen::MatrixXd& dft_orbitals);
  // Rebuilds ThreeCenterIntegrals, only works if the original basisobjects
  // still exist
  void Rebuild() { Fill(*_auxbasis, *_dftbasis, *_dft_orbitals); }

  void MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& AuxMatrix);

 private:
  // store vector of matrices
  std::vector<Eigen::MatrixXd> _matrix;

  // band summation indices
  int _mmin;
  int _mmax;
  int _nmin;
  int _nmax;
  int _ntotal;
  int _mtotal;
  int _basissize;

  const AOBasis* _auxbasis = nullptr;
  const AOBasis* _dftbasis = nullptr;
  const Eigen::MatrixXd* _dft_orbitals = nullptr;

  void FillBlock(std::vector<Eigen::MatrixXd>& matrix, const AOShell& auxshell,
                 const AOBasis& dftbasis, const Eigen::MatrixXd& dft_orbitals);
};

}  // namespace xtp
}  // namespace votca

#endif /* AOMATRIX_H */
