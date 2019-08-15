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
#ifndef _VOTCA_XTP_BSE_H
#define _VOTCA_XTP_BSE_H

#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {
struct BSE_Population;
template <int cqp, int cx, int cd, int cd2>
class BSE_OPERATOR;
typedef BSE_OPERATOR<1, 2, 1, 0> SingletOperator_TDA;
typedef BSE_OPERATOR<1, 0, 1, 0> TripletOperator_TDA;
template <class T>
class QMFragment;

class BSE {

 public:
  BSE(Logger& log, TCMatrix_gwbse& Mmn, const Eigen::MatrixXd& Hqp)
      : _log(log), _Mmn(Mmn), _Hqp(Hqp){};

  struct options {
    bool useTDA = true;
    int homo;
    int rpamin;
    int rpamax;
    int qpmin;
    int vmin;
    int cmax;
    int nmax = 5;             // number of eigenvectors to calculate
    bool davidson = true;     // use davidson to diagonalize the matrix
    bool matrixfree = false;  // use matrix free method
    std::string davidson_correction = "DPR";
    std::string davidson_ortho = "GS";
    std::string davidson_tolerance = "normal";
    std::string davidson_update = "safe";
    int davidson_maxiter = 50;
    double min_print_weight =
        0.5;  // minimium contribution for state to print it
  };

  void configure(const options& opt, const Eigen::VectorXd& DFTenergies);

  void Solve_singlets(Orbitals& orb) const;
  void Solve_triplets(Orbitals& orb) const;

  SingletOperator_TDA getSingletOperator_TDA() const;
  TripletOperator_TDA getTripletOperator_TDA() const;

  void Analyze_singlets(std::vector<QMFragment<BSE_Population> >& singlets,
                        const Orbitals& orb) const;
  void Analyze_triplets(std::vector<QMFragment<BSE_Population> >& triplets,
                        const Orbitals& orb) const;

 private:
  options _opt;

  struct Interaction {
    Eigen::VectorXd exchange_contrib;
    Eigen::VectorXd direct_contrib;
    Eigen::VectorXd qp_contrib;
  };

  Logger& _log;
  int _bse_vmax;
  int _bse_cmin;
  int _bse_size;
  int _bse_vtotal;
  int _bse_ctotal;

  Eigen::VectorXd _epsilon_0_inv;

  TCMatrix_gwbse& _Mmn;
  const Eigen::MatrixXd& _Hqp;

  tools::EigenSystem Solve_singlets_TDA() const;
  tools::EigenSystem Solve_singlets_BTDA() const;

  tools::EigenSystem Solve_triplets_TDA() const;
  tools::EigenSystem Solve_triplets_BTDA() const;

  void PrintWeights(const Eigen::VectorXd& weights) const;

  template <typename BSE_OPERATOR>
  void configureBSEOperator(BSE_OPERATOR& H) const;

  template <typename BSE_OPERATOR>
  tools::EigenSystem solve_hermitian(BSE_OPERATOR& H) const;

  template <typename BSE_OPERATOR_ApB, typename BSE_OPERATOR_AmB>
  tools::EigenSystem Solve_nonhermitian(BSE_OPERATOR_ApB& apb,
                                        BSE_OPERATOR_AmB&) const;

  void printFragInfo(const std::vector<QMFragment<BSE_Population> >& frags,
                     int state) const;
  void printWeights(int i_bse, double weight) const;
  void SetupDirectInteractionOperator(const Eigen::VectorXd& DFTenergies);

  Interaction Analyze_eh_interaction(const QMStateType& type,
                                     const Orbitals& orb) const;
  template <typename BSE_OPERATOR>
  Eigen::VectorXd Analyze_IndividualContribution(const QMStateType& type,
                                                 const Orbitals& orb,
                                                 const BSE_OPERATOR& H) const;
};
}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_BSE_H */
