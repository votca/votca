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

#ifndef VOTCA_XTP_ERIS_H
#define VOTCA_XTP_ERIS_H

#include <votca/xtp/fourcenter.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

/**
 * \brief Takes a density matrix and and an auxillary basis set and calculates
 * the electron repulsion integrals.
 *
 *
 *
 */
class ERIs {

 public:
  void Initialize(const AOBasis& dftbasis, const AOBasis& auxbasis);
  void Initialize_4c_small_molecule(const AOBasis& dftbasis);
  void Initialize_4c_screening(const AOBasis& dftbasis,
                               double eps);  // Pre-screening

  Mat_p_Energy CalculateERIs(const Eigen::MatrixXd& DMAT) const;
  Mat_p_Energy CalculateEXX(const Eigen::MatrixXd& DMAT) const;
  Mat_p_Energy CalculateEXX(const Eigen::MatrixXd& occMos,
                            const Eigen::MatrixXd& DMAT) const;
  Mat_p_Energy CalculateERIs_4c_small_molecule(
      const Eigen::MatrixXd& DMAT) const;
  Mat_p_Energy CalculateEXX_4c_small_molecule(
      const Eigen::MatrixXd& DMAT) const;

  Mat_p_Energy CalculateERIs_4c_direct(const AOBasis& dftbasis,
                                       const Eigen::MatrixXd& DMAT) const;

  int Removedfunctions() const { return _threecenter.Removedfunctions(); }

 private:
  bool _with_screening = false;
  double _screening_eps;
  Eigen::MatrixXd _diagonals;  // Square matrix containing <ab|ab> for all basis
                               // functions a, b

  void CalculateERIsDiagonals(const AOBasis& dftbasis);

  bool CheckScreen(double eps, const AOShell& shell_1, const AOShell& shell_2,
                   const AOShell& shell_3, const AOShell& shell_4) const;

  TCMatrix_dft _threecenter;
  FCMatrix _fourcenter;

  double CalculateEnergy(const Eigen::MatrixXd& DMAT,
                         const Eigen::MatrixXd& matrix_operator) const;
  template <bool transposed_block>
  void FillERIsBlock(Eigen::MatrixXd& ERIsCur, const Eigen::MatrixXd& DMAT,
                     const tensor4d& block, const AOShell& shell_1,
                     const AOShell& shell_2, const AOShell& shell_3,
                     const AOShell& shell_4) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ERIS_H
