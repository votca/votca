/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
#ifndef VOTCA_XTP_RPA_H
#define VOTCA_XTP_RPA_H

// Standard includes
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "logger.h"
#include "out_of_window_shift.h"
#include "screening_kernels.h"

namespace votca {
namespace xtp {
class TCMatrix_gwbse;

class RPA {
 public:
  RPA(Logger& log, const TCMatrix_gwbse& Mmn) : log_(log), Mmn_(Mmn){};

  void configure(Index homo, Index rpamin, Index rpamax) {
    homo_ = homo;
    rpamin_ = rpamin;
    rpamax_ = rpamax;
  }

  double getEta() const { return eta_; }

  /// How levels outside the GW window are shifted by UpdateRPAInputEnergies.
  void setOutOfWindowShift(const OutOfWindowShift& shift) {
    out_of_window_ = shift;
  }

  /// Times the dielectric-matrix builds in these timings (optional).
  void setTimings(DFTTimings* timings) { timings_ = timings; }

  Eigen::MatrixXd calculate_epsilon_i(double frequency) const {
    return calculate_epsilon<true>(frequency);
  }

  Eigen::MatrixXd calculate_epsilon_r(double frequency) const {
    return calculate_epsilon<false>(frequency);
  }

  Eigen::MatrixXd calculate_epsilon_r(std::complex<double> frequency) const;

  const Eigen::VectorXd& getRPAInputEnergies() const { return energies_; }

  void setRPAInputEnergies(const Eigen::VectorXd& rpaenergies) {
    energies_ = rpaenergies;
  }

  // calculates full RPA vector of energies from gwa and dftenergies and qpmin
  // RPA energies have three parts, lower than qpmin: dftenergies,between qpmin
  // and qpmax:gwa_energies,above:dftenergies+homo-lumo shift
  void UpdateRPAInputEnergies(const Eigen::VectorXd& dftenergies,
                              const Eigen::VectorXd& gwaenergies, Index qpmin);

  struct rpa_eigensolution {
    Eigen::VectorXd omega;    // Eigenvalues
    Eigen::MatrixXd XpY;      // Eigenvector components (X + Y)
    double ERPA_correlation;  // total correlation energy
  };

  rpa_eigensolution Diagonalize_H2p() const;

 private:
  Index homo_;  // HOMO index with respect to dft energies
  Index rpamin_;
  Index rpamax_;
  const double eta_ = 0.0001;

  Eigen::VectorXd energies_;

  Logger& log_;
  const TCMatrix_gwbse& Mmn_;
  DFTTimings* timings_ = nullptr;

  // Calls use() with the virtual rows of hole slice m_level.
  void VisitHoleVirtualRows(Index m_level,
                            const WeightedGram::RowsVisitor& use) const;
  // sum_v M_v^T diag(w_v) M_v over the hole slices v, M_v their virtual rows
  Eigen::MatrixXd ResponseSum(const WeightedGram::WeightsFn& weights) const;

  template <bool imag>
  Eigen::MatrixXd calculate_epsilon(double frequency) const;

  Eigen::VectorXd Calculate_H2p_AmB() const;
  Eigen::MatrixXd Calculate_H2p_ApB() const;
  SymmetricEigenSystem Diagonalize_H2p_C(const Eigen::MatrixXd& C) const;

  void ShiftUncorrectedEnergies(const Eigen::VectorXd& dftenergies, Index qpmin,
                                Index gwsize);

  OutOfWindowShift out_of_window_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_RPA_H
