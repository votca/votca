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
#ifndef _VOTCA_XTP_RPA_H
#define _VOTCA_XTP_RPA_H
#include "votca/xtp/logger.h"
#include <vector>
#include <votca/tools/eigensystem.h>
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {
class TCMatrix_gwbse;

class RPA {
 public:
  RPA(Logger& log, const TCMatrix_gwbse& Mmn) : _log(log), _Mmn(Mmn){};

  void configure(Index homo, Index rpamin, Index rpamax) {
    _homo = homo;
    _rpamin = rpamin;
    _rpamax = rpamax;
  }

  Eigen::MatrixXd calculate_epsilon_i(double frequency) const {
    return calculate_epsilon<true>(frequency);
  }

  Eigen::MatrixXd calculate_epsilon_r(double frequency) const {
    return calculate_epsilon<false>(frequency);
  }

  const Eigen::VectorXd& getRPAInputEnergies() const { return _energies; }

  void setRPAInputEnergies(const Eigen::VectorXd& rpaenergies) {
    _energies = rpaenergies;
  }

  // calculates full RPA vector of energies from gwa and dftenergies and qpmin
  // RPA energies have three parts, lower than qpmin: dftenergies,between qpmin
  // and qpmax:gwa_energies,above:dftenergies+homo-lumo shift
  void UpdateRPAInputEnergies(const Eigen::VectorXd& dftenergies,
                              const Eigen::VectorXd& gwaenergies, Index qpmin);

  struct rpa_eigensolution {
    Eigen::VectorXd _omega;  // Eigenvalues
    Eigen::MatrixXd _XpY;    // Eigenvector components (X + Y)
  };

  rpa_eigensolution Diagonalize_H2p() const;

 private:
  Index _homo;  // HOMO index with respect to dft energies
  Index _rpamin;
  Index _rpamax;
  const double _eta = 0.0001;

  Eigen::VectorXd _energies;

  Logger& _log;
  const TCMatrix_gwbse& _Mmn;

  template <bool imag>
  Eigen::MatrixXd calculate_epsilon(double frequency) const;

  Eigen::VectorXd Calculate_H2p_AmB() const;
  Eigen::MatrixXd Calculate_H2p_ApB() const;
  tools::EigenSystem Diagonalize_H2p_C(const Eigen::MatrixXd& C) const;
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_RPA_RPA_H */
