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

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <votca/tools/constants.h>
#include <votca/xtp/sigma_base.h>

#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

Eigen::MatrixXd Sigma_base::CalcExchange() const {

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
  int gwsize = _Mmn.auxsize();
  int occlevel = _opt.homo - _opt.rpamin + 1;
  int qpmin = _opt.qpmin - _opt.rpamin;
#pragma omp parallel for schedule(dynamic)
  for (int gw_level1 = 0; gw_level1 < _qptotal; gw_level1++) {
    const Eigen::MatrixXd& Mmn1 = _Mmn[gw_level1 + qpmin];
    for (int gw_level2 = gw_level1; gw_level2 < _qptotal; gw_level2++) {
      const Eigen::MatrixXd& Mmn2 = _Mmn[gw_level2 + qpmin];
      double sigma_x = -(Mmn1.block(0, 0, occlevel, gwsize)
                             .cwiseProduct(Mmn2.block(0, 0, occlevel, gwsize)))
                            .sum();
      result(gw_level1, gw_level2) = sigma_x;
      result(gw_level2, gw_level1) = sigma_x;
    }
  }
  return result;
}

}  // namespace xtp
}  // namespace votca
