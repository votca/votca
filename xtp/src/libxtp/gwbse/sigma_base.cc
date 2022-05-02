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

// Standard includes
#include <cmath>

// Third party includes
#include <boost/math/constants/constants.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/sigma_base.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

Eigen::MatrixXd Sigma_base::CalcExchangeMatrix() const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
  Index occlevel = opt_.homo - opt_.rpamin + 1;
  Index qpmin = opt_.qpmin - opt_.rpamin;
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level1 = 0; gw_level1 < qptotal_; gw_level1++) {
    const Eigen::MatrixXd& Mmn1 = Mmn_[gw_level1 + qpmin];
    for (Index gw_level2 = gw_level1; gw_level2 < qptotal_; gw_level2++) {
      const Eigen::MatrixXd& Mmn2 = Mmn_[gw_level2 + qpmin];
      double sigma_x =
          -(Mmn1.topRows(occlevel).cwiseProduct(Mmn2.topRows(occlevel))).sum();
      result(gw_level2, gw_level1) = sigma_x;
    }
  }
  result = result.selfadjointView<Eigen::Lower>();
  return result;
}

Eigen::VectorXd Sigma_base::CalcCorrelationDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(qptotal_);
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; gw_level++) {
    result(gw_level) =
        CalcCorrelationDiagElement(gw_level, frequencies[gw_level]);
  }
  return result;
}

Eigen::MatrixXd Sigma_base::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level1 = 0; gw_level1 < qptotal_; gw_level1++) {
    for (Index gw_level2 = gw_level1 + 1; gw_level2 < qptotal_; gw_level2++) {
      double sigma_c = CalcCorrelationOffDiagElement(
          gw_level1, gw_level2, frequencies[gw_level1], frequencies[gw_level2]);
      result(gw_level2, gw_level1) = sigma_c;
    }
  }
  result = result.selfadjointView<Eigen::Lower>();
  return result;
}

}  // namespace xtp
}  // namespace votca
