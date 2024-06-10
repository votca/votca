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
#ifndef VOTCA_XTP_BSE_H
#define VOTCA_XTP_BSE_H

// Local VOTCA includes
#include "logger.h"
#include "orbitals.h"
#include "qmstate.h"
#include "threecenter.h"

namespace votca {
namespace xtp {
struct BSE_Population;
template <Index cqp, Index cx, Index cd, Index cd2>
class BSE_OPERATOR;
typedef BSE_OPERATOR<1, 2, 1, 0> SingletOperator_TDA;
typedef BSE_OPERATOR<1, 0, 1, 0> TripletOperator_TDA;
template <class T>
class QMFragment;

class BSE {

 public:
  //  BSE(Logger& log, TCMatrix_gwbse& Mmn, const Eigen::MatrixXd& Hqp_in)
  //    :  log_(log),  Mmn_(Mmn),  Hqp_in_(Hqp_in){};
  BSE(Logger& log, TCMatrix_gwbse& Mmn) : log_(log), Mmn_(Mmn) {};

  struct options {
    bool useTDA;
    Index homo;
    Index rpamin;
    Index rpamax;
    Index qpmin;
    Index qpmax;
    Index vmin;
    Index cmax;
    Index nmax;  // number of eigenvectors to calculat
    std::string davidson_correction;
    std::string davidson_tolerance;
    std::string davidson_update;
    Index davidson_maxiter;
    double min_print_weight;  // minimium contribution for state to print it
    bool use_Hqp_offdiag;
    Index max_dyn_iter;
    double dyn_tolerance;
  };

  void configure(const options& opt, const Eigen::VectorXd& RPAEnergies,
                 const Eigen::MatrixXd& Hqp_in);

  void Solve_singlets(Orbitals& orb) const;
  void Solve_triplets(Orbitals& orb) const;

  Eigen::MatrixXd getHqp() const { return Hqp_; };

  SingletOperator_TDA getSingletOperator_TDA() const;
  TripletOperator_TDA getTripletOperator_TDA() const;

  void Analyze_singlets(std::vector<QMFragment<BSE_Population> > fragments,
                        const Orbitals& orb) const;
  void Analyze_triplets(std::vector<QMFragment<BSE_Population> > fragments,
                        const Orbitals& orb) const;

  void Perturbative_DynamicalScreening(const QMStateType& type, Orbitals& orb);

 private:
  options opt_;

  struct Interaction {
    Eigen::VectorXd exchange_contrib;
    Eigen::VectorXd direct_contrib;
    Eigen::VectorXd qp_contrib;
  };

  struct ExpectationValues {
    Eigen::VectorXd direct_term;
    Eigen::VectorXd cross_term;
  };

  Logger& log_;
  Index bse_vmax_;
  Index bse_cmin_;
  Index bse_size_;
  Index bse_vtotal_;
  Index bse_ctotal_;

  Index max_dyn_iter_;
  double dyn_tolerance_;

  Eigen::VectorXd epsilon_0_inv_;

  TCMatrix_gwbse& Mmn_;
  Eigen::MatrixXd Hqp_;

  tools::EigenSystem Solve_singlets_TDA() const;
  tools::EigenSystem Solve_singlets_BTDA() const;

  tools::EigenSystem Solve_triplets_TDA() const;
  tools::EigenSystem Solve_triplets_BTDA() const;

  void PrintWeights(const Eigen::VectorXd& weights) const;

  template <typename BSE_OPERATOR>
  void configureBSEOperator(BSE_OPERATOR& H) const;

  template <typename BSE_OPERATOR>
  tools::EigenSystem solve_hermitian(BSE_OPERATOR& h) const;

  template <typename BSE_OPERATOR_ApB, typename BSE_OPERATOR_AmB>
  tools::EigenSystem Solve_nonhermitian(BSE_OPERATOR_ApB& apb,
                                        BSE_OPERATOR_AmB&) const;

  template <typename BSE_OPERATOR_A, typename BSE_OPERATOR_B>
  tools::EigenSystem Solve_nonhermitian_Davidson(BSE_OPERATOR_A& Aop,
                                                 BSE_OPERATOR_B& Bop) const;

  void printFragInfo(const std::vector<QMFragment<BSE_Population> >& frags,
                     Index state) const;
  void printWeights(Index i_bse, double weight) const;
  void SetupDirectInteractionOperator(const Eigen::VectorXd& DFTenergies,
                                      double energy);

  Eigen::MatrixXd AdjustHqpSize(const Eigen::MatrixXd& Hqp_in,
                                const Eigen::VectorXd& RPAInputEnergies);

  Interaction Analyze_eh_interaction(const QMStateType& type,
                                     const Orbitals& orb) const;
  template <typename BSE_OPERATOR>
  ExpectationValues ExpectationValue_Operator(const QMStateType& type,
                                              const Orbitals& orb,
                                              const BSE_OPERATOR& H) const;

  template <typename BSE_OPERATOR>
  ExpectationValues ExpectationValue_Operator_State(
      const QMState& state, const Orbitals& orb, const BSE_OPERATOR& H) const;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSE_H
