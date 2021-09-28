
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
#ifndef VOTCA_XTP_BSECOUPLING_H
#define VOTCA_XTP_BSECOUPLING_H

// Local VOTCA includes
#include "couplingbase.h"
#include "qmstate.h"

namespace votca {
namespace xtp {

/**
 * \brief Evaluates electronic coupling elements
 *
 * J. Wehner,B. Baumeier,
 * JCTC DOI: 10.1021/acs.jctc.6b00935
 *
 */

class BSECoupling : public CouplingBase {
 public:
  void Initialize(tools::Property& options) override;
  std::string Identify() const { return "bsecoupling"; }

  void Addoutput(tools::Property& type_summary, const Orbitals& orbitalsA,
                 const Orbitals& orbitalsB) const override;

  /**
   * \brief evaluates electronic couplings
   *
   * @param  orbitalsA molecular orbitals of molecule A
   * @param  orbitalsB molecular orbitals of molecule B
   * @param  orbitalsAB molecular orbitals of the dimer AB
   */
  void CalculateCouplings(const Orbitals& orbitalsA, const Orbitals& orbitalsB,
                          const Orbitals& orbitalsAB) override;

 private:
  void WriteToProperty(tools::Property& summary, const QMState& stateA,
                       const QMState& stateB) const;

  double getSingletCouplingElement(Index levelA, Index levelB,
                                   Index methodindex) const;

  double getTripletCouplingElement(Index levelA, Index levelB,
                                   Index methodindex) const;

  Eigen::MatrixXd SetupCTStates(Index bseA_vtotal, Index bseB_vtotal,
                                Index bseAB_vtotal, Index bseAB_ctotal,
                                const Eigen::MatrixXd& A_AB,
                                const Eigen::MatrixXd& B_AB) const;

  Eigen::MatrixXd ProjectFrenkelExcitons(const Eigen::MatrixXd& BSE_Coeffs,
                                         const Eigen::MatrixXd& X_AB,
                                         Index bseX_vtotal, Index bseX_ctotal,
                                         Index bseAB_vtotal,
                                         Index bseAB_ctotal) const;

  template <class BSE_OPERATOR>
  std::array<Eigen::MatrixXd, 2> ProjectExcitons(Eigen::MatrixXd& FE_AB,
                                                 Eigen::MatrixXd& CTStates,
                                                 BSE_OPERATOR H) const;
  template <class BSE_OPERATOR>
  Eigen::MatrixXd CalcJ_dimer(BSE_OPERATOR& H,
                              Eigen::MatrixXd& projection) const;

  Eigen::MatrixXd OrthogonalizeCTs(Eigen::MatrixXd& FE_AB,
                                   Eigen::MatrixXd& CTStates) const;

  Eigen::MatrixXd Fulldiag(const Eigen::MatrixXd& J_dimer) const;

  Eigen::MatrixXd Perturbation(const Eigen::MatrixXd& J_dimer) const;

  std::array<Eigen::MatrixXd, 2> JAB_singlet;
  std::array<Eigen::MatrixXd, 2> JAB_triplet;

  bool doTriplets_;
  bool doSinglets_;
  bool output_perturbation_;
  Index levA_;
  Index levB_;
  Index occA_;
  Index unoccA_;
  Index occB_;
  Index unoccB_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSECOUPLING_H
