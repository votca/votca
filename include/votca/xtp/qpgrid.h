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

 public:
  QPGrid(Logger& log,
         Sigma_base& sigma,
         const Eigen::MatrixXd& vxc_diag,
         const Eigen::VectorXd& dft_energies,
         const Eigen::VectorXd& sigma_x_diag)
      : _log(log),
        _sigma(sigma),
        _vxc_diag(vxc_diag),
        _dft_energies(dft_energies),
        _sigma_x_diag(sigma_x_diag) {}

  struct options {
    Index homo = 0;
    Index qpmin = 0;
    Index qpmax = 0;
    Index steps = 201;
    double range = 1.0;
  };

  void configure(options opt) {
    _opt = opt;
    _qptotal = opt.qpmax - opt.qpmin + 1;
    _steps = opt.steps;
    _range = opt.range;
    _delta = 2.0 * _range / (_steps - 1);
  }

  Eigen::VectorXd Evaluate(Eigen::VectorXd frequencies);

 private:
  Logger& _log;
  options _opt;
  Index _qptotal = 0;
  Index _steps = 0; // Number of grid points
  double _range = 0.0; // Distance left/right from the origin
  double _delta = 0.0; // Step size
  Sigma_base& _sigma;

  const Eigen::VectorXd& _vxc_diag;
  const Eigen::VectorXd& _dft_energies;
  const Eigen::VectorXd& _sigma_x_diag;

};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_QPGRID_H */
