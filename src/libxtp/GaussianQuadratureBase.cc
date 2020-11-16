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

// Local VOTCA includes
#include "votca/xtp/GaussianQuadratureBase.h"

namespace votca {
namespace xtp {
void GaussianQuadratureBase::CheckOrder(
    Index order, const std::map<Index, Eigen::VectorXd>& map) const {
  if (map.count(order) == 0) {
    std::string keys = "{ ";
    for (const auto& pair : map) {
      keys += std::to_string(pair.first) + " ";
    }
    keys += "}";
    throw std::invalid_argument("Order " + std::to_string(order) + " not in " +
                                keys +
                                ". Please "
                                " select one of these numbers");
  }
}

}  // namespace xtp
}  // namespace votca
