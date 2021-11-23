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
#ifndef VOTCA_XTP_MULTIPOLE_H
#define VOTCA_XTP_MULTIPOLE_H

#include "eigen.h"
#include "eigen3/Eigen/Dense"

namespace votca {
namespace xtp {

struct Multipole {
  Index rank = 0;
  double charge = 0.0;
  Eigen::Vector3d dipole = Eigen::Vector3d::Zero();
  Eigen::Matrix3d quadrupole = Eigen::Matrix3d::Zero();
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EANALYZE_H
