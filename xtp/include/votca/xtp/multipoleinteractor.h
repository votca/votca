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

#pragma once
#ifndef VOTCA_XTP_MULTIPOLEINTERACTOR_H
#define VOTCA_XTP_MULTIPOLEINTERACTOR_H

#include <array>
#include <boost/math/constants/constants.hpp>
#include <vector>

#include "multipole.h"

namespace votca {
namespace xtp {

class MultipoleInteractor {
 public:
  MultipoleInteractor(double alpha, double thole_damping);

  std::array<double, 5> orientationDependence(const Multipole& mp1,
                                              const Multipole& mp2,
                                              const Eigen::Vector3d& dr) const;

  std::array<Eigen::Vector3d, 3> orientationDependence(
      const Multipole& mp1, const Eigen::Vector3d& dr) const;

  std::array<double, 5> erfDistanceDependence(const Eigen::vector3d& dr) const;
  std::array<double, 5> erfcDistanceDependence(const Eigen::vector3d& dr) const;

  std::array<double, 5> tholeDamping(const Eigen::Matrix3d& pol1,
                                     const Eigen::Matrix3d& pol2,
                                     const Eigen::Vector3d& dr) const;

 private:
  constexpr double pi = boost::math::constants::pi<double>();
  constexpr double rSqrtPi = 1.0 / std::sqrt(pi);

  double a1, a2, a3, a4, a5;     // alpha (splitting param) and its powers
  double thole, thole2, thole3;  // thole damping params and its powers
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EANALYZE_H
