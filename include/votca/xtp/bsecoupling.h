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

#include <votca/xtp/couplingbase.h>
#include <votca/xtp/qmstate.h>

#pragma once
#ifndef VOTCA_XTP_BSECOUPLING_H
#define VOTCA_XTP_BSECOUPLING_H

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
  void Initialize(tools::Property& options);
  std::string Identify() const { return "bsecoupling"; }

  Eigen::MatrixXd getJAB_singletstorage() const {
    return (_output_perturbation ? JAB_singlet[0] : JAB_singlet[1]);
  }

  Eigen::MatrixXd getJAB_tripletstorage() const {
    return (_output_perturbation ? JAB_triplet[0] : JAB_triplet[1]);
  }
  void Addoutput(tools::Property& type_summary, const Orbitals& orbitalsA,
                 const Orbitals& orbitalsB) const;

  /**
   * \brief evaluates electronic couplings
   *
   * @param _orbitalsA molecular orbitals of molecule A
   * @param _orbitalsB molecular orbitals of molecule B
   * @param _orbitalsAB molecular orbitals of the dimer AB
   */
  void CalculateCouplings(const Orbitals& orbitalsA, const Orbitals& orbitalsB,
                          const Orbitals& orbitalsAB);

 private:
  void WriteToProperty(tools::Property& summary, const QMState& stateA,
                       const QMState& stateB) const;

  double getSingletCouplingElement(int levelA, int levelB,
                                   int methodindex) const;

  double getTripletCouplingElement(int levelA, int levelB,
                                   int methodindex) const;

  Eigen::MatrixXd SetupCTStates(int bseA_vtotal, int bseB_vtotal,
                                int bseAB_vtotal, int bseAB_ctotal,
                                const Eigen::MatrixXd& A_AB,
                                const Eigen::MatrixXd& B_AB) const;

  Eigen::MatrixXd ProjectFrenkelExcitons(const Eigen::MatrixXd& BSE_Coeffs,
                                         const Eigen::MatrixXd& X_AB,
                                         int bseX_vtotal, int bseX_ctotal,
                                         int bseAB_vtotal,
                                         int bseAB_ctotal) const;

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

  bool _doTriplets = false;
  bool _doSinglets = false;
  bool _output_perturbation = true;
  int _levA = std::numeric_limits<int>::max();
  int _levB = std::numeric_limits<int>::max();
  int _occA = 5;
  int _unoccA = 5;
  int _occB = 5;
  int _unoccB = 5;

  // Concurrent jobs in the GPU
  int _max_gpu_streams;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSECOUPLING_H
