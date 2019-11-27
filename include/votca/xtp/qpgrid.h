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
  QPGrid(Sigma_base& sigma) : _sigma(sigma) {}

  struct options {
    Index homo = 0;
    Index qpmin = 0;
    Index qpmax = 0;
    Index steps = 201;
    double spacing = 0.01;
  };

  void configure(options opt) { _opt = opt; }

  void setQPOffset(const Eigen::VectorXd& qpOffset) {
    _qpOffset = qpOffset;
  }

  Eigen::VectorXd Evaluate(const Eigen::VectorXd frequencies);

 private:
  Sigma_base& _sigma;
  options _opt;
  Eigen::VectorXd _qpOffset;
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_QPGRID_H */
