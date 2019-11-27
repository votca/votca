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

  // This class solves the QP equation for the variable QP energy E_QP:
  //   \Sigma_c(E_QP) = E_QP - E_intercept.
  // Here, E_intercept denotes the intercept of the line function on the
  // right-hand-side of the equation:
  //   E_intercept = E_KS + \Sigma_x - v_xc,
  // with
  //   E_KS the Kohn-Sham energy,
  //   \Sigma_x the exchange part of the self-energy,
  //   v_xc the exchange-correlation potential.

 public:
  QPGrid(Sigma_base& sigma) : _sigma(sigma) {}

  struct options {
    Index homo = 0;
    Index qpmin = 0;
    Index qpmax = 0;
    Index steps = 201;  // Number of grid points
    double spacing = 0.01;  // Spacing of grid points in Ha
  };

  void configure(options opt) { _opt = opt; }

  void setEnergyIntercept(const Eigen::VectorXd& intercept) {
    _energy_intercept = intercept;
  }

  Eigen::VectorXd Evaluate(const Eigen::VectorXd frequencies);

 private:
  Sigma_base& _sigma;
  options _opt;
  Eigen::VectorXd _energy_intercept;
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_QPGRID_H */
