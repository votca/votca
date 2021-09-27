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
#ifndef VOTCA_XTP_ERIS_H
#define VOTCA_XTP_ERIS_H

// Local VOTCA includes
#include "threecenter.h"

namespace votca {
namespace xtp {

/**
 * \brief Takes a density matrix and and an auxiliary basis set and calculates
 * the electron repulsion integrals.
 *
 */
class ERIs {

 public:
  void Initialize(const AOBasis& dftbasis, const AOBasis& auxbasis);
  void Initialize_4c(const AOBasis& dftbasis);

  Eigen::MatrixXd CalculateERIs_3c(const Eigen::MatrixXd& DMAT) const;

  std::array<Eigen::MatrixXd, 2> CalculateERIs_EXX_3c(
      const Eigen::MatrixXd& occMos, const Eigen::MatrixXd& DMAT) const {
    std::array<Eigen::MatrixXd, 2> result;
    result[0] = CalculateERIs_3c(DMAT);
    if (occMos.rows() > 0 && occMos.cols() > 0) {
      assert(occMos.rows() == DMAT.rows() && "occMos.rows()==DMAT.rows()");
      result[1] = CalculateEXX_mos(occMos);
    } else {
      result[1] = CalculateEXX_dmat(DMAT);
    }
    return result;
  }

  Eigen::MatrixXd CalculateERIs_4c(const Eigen::MatrixXd& DMAT,
                                   double error) const {
    return Compute4c<false>(DMAT, error)[0];
  }

  std::array<Eigen::MatrixXd, 2> CalculateERIs_EXX_4c(
      const Eigen::MatrixXd& DMAT, double error) const {
    return Compute4c<true>(DMAT, error);
  }

  Index Removedfunctions() const { return _threecenter.Removedfunctions(); }

  static double CalculateEnergy(const Eigen::MatrixXd& DMAT,
                                const Eigen::MatrixXd& matrix_operator) {
    return matrix_operator.cwiseProduct(DMAT).sum();
  }

 private:
  std::vector<libint2::Shell> basis_;
  std::vector<Index> starts_;

  std::vector<std::vector<Index>> shellpairs_;
  std::vector<std::vector<libint2::ShellPair>> shellpairdata_;
  Index maxnprim_;
  Index maxL_;

  Eigen::MatrixXd CalculateEXX_dmat(const Eigen::MatrixXd& DMAT) const;
  Eigen::MatrixXd CalculateEXX_mos(const Eigen::MatrixXd& occMos) const;

  std::vector<std::vector<libint2::ShellPair>> ComputeShellPairData(
      const std::vector<libint2::Shell>& basis,
      const std::vector<std::vector<Index>>& shellpairs) const;

  Eigen::MatrixXd ComputeSchwarzShells(const AOBasis& dftbasis) const;
  Eigen::MatrixXd ComputeShellBlockNorm(const Eigen::MatrixXd& dmat) const;

  template <bool with_exchange>
  std::array<Eigen::MatrixXd, 2> Compute4c(const Eigen::MatrixXd& dmat,
                                           double error) const;

  TCMatrix_dft _threecenter;

  Eigen::MatrixXd schwarzscreen_;  // Square matrix containing <ab|ab> for all
                                   // shells
};                                 // namespace xtp

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ERIS_H
