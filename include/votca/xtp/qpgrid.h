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
#ifndef _VOTCA_XTP_QPGRID_H
#define _VOTCA_XTP_QPGRID_H
#include "votca/xtp/eigen.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/sigma_base.h"

namespace votca {
namespace xtp {

class QPGrid {

  // This class solves the QP equation for the variable QP energy '\omega':
  //   \Sigma_c(\omega) = \omega - E_intercept.
  // Every solution '\omega' is associated to a pole in the spectral function
  // 'A(\omega)'. 'E_intercept' denotes the intercept of the line function on
  // the right-hand-side of the equation:
  //   E_intercept = E_KS + \Sigma_x - v_xc,
  // with
  //   'E_KS' the Kohn-Sham energy,
  //   '\Sigma_x' the exchange part of the self-energy,
  //   'v_xc' the exchange-correlation potential.
  //
  // The QP energy is defined as the solution '\omega*' for which the pole
  // weight 'Z' is maximized. Here, we assume that the pole with the highest
  // pole weight corresponds to the coherent part of the spectral function
  // 'A(/omega)'.

 public:
  QPGrid(Sigma_base& sigma) : _sigma(sigma) {}

  struct options {
    Index homo = 0;
    Index qpmin = 0;
    Index qpmax = 0;
    Index steps = 201;      // Number of grid points
    double spacing = 0.01;  // Spacing of grid points in Ha
  };

  void configure(options opt) { _opt = opt; }

  void setEnergyIntercept(const Eigen::VectorXd& intercept) {
    _energy_intercept = intercept;
  }

  Eigen::VectorXd Evaluate(const Eigen::VectorXd& gwa_energies);

 private:
  Sigma_base& _sigma;
  options _opt;
  Eigen::VectorXd _energy_intercept;

  // 'FindQPEnergy' finds the best estimate of the QP energy 'qp_energy' on the
  // grid with grid points 'freq'. The self-energy correlation evaluated at the
  // grid points is contained by 'sigc'. Returns 'true' if a solution to the QP
  // equation was found and 'false' otherwise.
  //
  // The solutions to the QP equation are found by solving the root-finding
  // problem associated to the target function 'f_targ':
  //   f_targ(\omega) = \Sigma_c(\omega) + E_intercept - \omega
  // The target function evaluated at the grid points is contained by 'targ'.
  bool FindQPEnergy(const Eigen::VectorXd& freq, const Eigen::VectorXd& sigc,
                    const Eigen::VectorXd& targ, double& qp_energy) const;
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_QPGRID_H */
