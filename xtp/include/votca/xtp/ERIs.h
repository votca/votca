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
  /// pair_threshold: shell pairs whose three-centre integrals are bounded
  /// by less than this are not stored (see TCMatrix_dft); 0 keeps all.
  void Initialize(const AOBasis& dftbasis, const AOBasis& auxbasis,
                  double pair_threshold = 0.0);
  void Initialize_4c(const AOBasis& dftbasis);

  Eigen::MatrixXd CalculateERIs_3c(const Eigen::MatrixXd& DMAT) const;

  std::array<Eigen::MatrixXd, 2> CalculateERIs_EXX_3c(
      const Eigen::MatrixXd& occMos, const Eigen::MatrixXd& DMAT) const {
    std::array<Eigen::MatrixXd, 2> result;
    result[0] = CalculateERIs_3c(DMAT);
    result[1] = CalculateEXX_3c(occMos, DMAT);
    return result;
  }

  /// Exchange matrix only: from the occupied MOs if given, otherwise from
  /// the density matrix.
  Eigen::MatrixXd CalculateEXX_3c(const Eigen::MatrixXd& occMos,
                                  const Eigen::MatrixXd& DMAT) const {
    if (occMos.rows() > 0 && occMos.cols() > 0) {
      assert(occMos.rows() == DMAT.rows() && "occMos.rows()==DMAT.rows()");
      return CalculateEXX_mos(occMos);
    }
    return CalculateEXX_dmat(DMAT);
  }

  /// Number of auxiliary functions the 3c tensor is stored for.
  Index AuxSize() const { return threecenter_.size(); }
  /// Stored basis-function pairs of the 3c tensor and all N(N+1)/2 pairs.
  Index StoredPairs() const { return threecenter_.StoredPairs(); }
  Index AllPairs() const { return threecenter_.AllPairs(); }
  /// Wall time spent on the aux Coulomb metric and its inverse square root.
  double MetricSeconds() const { return threecenter_.MetricSeconds(); }

  Eigen::MatrixXd CalculateERIs_4c(const Eigen::MatrixXd& DMAT,
                                   double error) const {
    return Compute4c<false>(DMAT, error)[0];
  }

  std::array<Eigen::MatrixXd, 2> CalculateERIs_EXX_4c(
      const Eigen::MatrixXd& DMAT, double error) const {
    return Compute4c<true>(DMAT, error);
  }

  Index Removedfunctions() const { return threecenter_.Removedfunctions(); }

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

  /// Relative size below which eigenvalues of a density matrix are treated
  /// as zero when it is factorised for the exchange build.
  static constexpr double kDensityRankCutoff = 1e-12;
  Eigen::MatrixXd CalculateEXX_dmat(const Eigen::MatrixXd& DMAT) const;
  Eigen::MatrixXd CalculateEXX_mos(const Eigen::MatrixXd& occMos) const;
  /// -sum_P (B_P X)(B_P X)^T + sum_P (B_P Y)(B_P Y)^T for factors = [X Y],
  /// X being the first npos columns.
  Eigen::MatrixXd ExchangeFromFactors(const Eigen::MatrixXd& factors,
                                      Index npos) const;
  /// Rows of B_P X collected per thread before one rank update.
  static constexpr Index kExchangeBatchRows = 1024;

  std::vector<std::vector<libint2::ShellPair>> ComputeShellPairData(
      const std::vector<libint2::Shell>& basis,
      const std::vector<std::vector<Index>>& shellpairs) const;

  Eigen::MatrixXd ComputeSchwarzShells(const AOBasis& dftbasis) const;
  Eigen::MatrixXd ComputeShellBlockNorm(const Eigen::MatrixXd& dmat) const;

  template <bool with_exchange>
  std::array<Eigen::MatrixXd, 2> Compute4c(const Eigen::MatrixXd& dmat,
                                           double error) const;

  TCMatrix_dft threecenter_;

  Eigen::MatrixXd schwarzscreen_;  // Square matrix containing <ab|ab> for all
                                   // shells
};  // namespace xtp

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ERIS_H
