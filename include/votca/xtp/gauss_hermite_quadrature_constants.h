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

#ifndef __VOTCA_XTP_GAUS_HERMITE_QUADRATURE_CONSTANTS_H
#define __VOTCA_XTP_GAUS_HERMITE_QUADRATURE_CONSTANTS_H

#include "votca/xtp/eigen.h"
#include <map>
#include <stdexcept>
#include <string>

namespace votca {
namespace xtp {

class Gauss_Hermite_Quadrature_Constants {
 public:
  const Eigen::VectorXd &getPoints(Index order);

  const Eigen::VectorXd &getAdaptedWeights(Index order);

 private:
  bool _filled_Points = false;
  bool _filled_AdaptedWeights = false;

  std::map<Index, Eigen::VectorXd> _map_points;
  std::map<Index, Eigen::VectorXd> _map_AdaptedWeights;

  void FillPoints();
  void FillAdaptedWeights();
};
}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_GAUS_HERMITE_QUADRATURE_CONSTANTS_H */
