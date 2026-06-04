/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#ifndef VOTCA_XTP_BSE_UKS_H
#define VOTCA_XTP_BSE_UKS_H

#include "logger.h"
#include "orbitals.h"
#include "qmfragment.h"
#include "threecenter.h"

namespace votca {
namespace xtp {

struct BSE_Population;

class BSE_UKS {
 public:
  struct options {
    bool useTDA;
    Index rpamin;
    Index rpamax;
    Index qpmin;
    Index qpmax;
    Index vmin;
    Index cmax;
    Index nmax;
    std::string davidson_correction;
    std::string davidson_tolerance;
    std::string davidson_update;
    Index davidson_maxiter;
    double min_print_weight;
    bool use_Hqp_offdiag;
    Index max_dyn_iter;
    double dyn_tolerance;
  };

  BSE_UKS(Logger& log, const TCMatrix_gwbse_spin& Mmn)
      : log_(log), Mmn_raw_(Mmn), Mmn_(Mmn) {}

  void configure_with_precomputed_screening(
      const options& opt, Index homo_alpha, Index homo_beta,
      const Eigen::VectorXd& RPAInputEnergiesAlpha,
      const Eigen::VectorXd& RPAInputEnergiesBeta,
      const Eigen::MatrixXd& Hqp_alpha_in, const Eigen::MatrixXd& Hqp_beta_in,
      const Eigen::VectorXd& epsilon_0_inv,
      const Eigen::MatrixXd& epsilon_eigenvectors);

  void Solve_excitons_uks(Orbitals& orb) const;

  void Analyze_excitons_uks(std::vector<QMFragment<BSE_Population>> fragments,
                            const Orbitals& orb) const;

  void Perturbative_DynamicalScreening(Orbitals& orb);

 private:
  struct ExpectationValues {
    Eigen::VectorXd direct_term;
    Eigen::VectorXd cross_term;
  };

  Logger& log_;
  const TCMatrix_gwbse_spin& Mmn_raw_;
  TCMatrix_gwbse_spin Mmn_;
  options opt_;

  Index homo_alpha_ = 0;
  Index homo_beta_ = 0;

  Index alpha_vtotal_ = 0;
  Index alpha_ctotal_ = 0;
  Index alpha_size_ = 0;
  Index beta_vtotal_ = 0;
  Index beta_ctotal_ = 0;
  Index beta_size_ = 0;

  Eigen::VectorXd epsilon_0_inv_;
  Eigen::MatrixXd Hqp_alpha_;
  Eigen::MatrixXd Hqp_beta_;

  Eigen::MatrixXd AdjustHqpSize(const Eigen::MatrixXd& Hqp_in,
                                const Eigen::VectorXd& RPAInputEnergies,
                                Index homo) const;

  void SetupDirectInteractionOperator(
      const Eigen::VectorXd& RPAInputEnergiesAlpha,
      const Eigen::VectorXd& RPAInputEnergiesBeta, double energy);

  template <typename BSE_OPERATOR>
  void configureBSEOperator(BSE_OPERATOR& H) const;

  template <typename BSE_OPERATOR>
  ExpectationValues ExpectationValue_Operator(const Orbitals& orb,
                                              const BSE_OPERATOR& H) const;

  template <typename BSE_OPERATOR>
  ExpectationValues ExpectationValue_Operator_State(
      Index state, const Orbitals& orb, const BSE_OPERATOR& H) const;

  tools::EigenSystem Solve_excitons_uks_TDA() const;
  tools::EigenSystem Solve_excitons_uks_BTDA() const;

  template <typename BSE_OPERATOR>
  tools::EigenSystem solve_hermitian(BSE_OPERATOR& h) const;

  template <typename BSE_OPERATOR_A, typename BSE_OPERATOR_B>
  tools::EigenSystem Solve_nonhermitian_Davidson(BSE_OPERATOR_A& Aop,
                                                 BSE_OPERATOR_B& Bop) const;

  void PrintWeightsUKS(const Eigen::VectorXd& coeffs) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSE_UKS_H