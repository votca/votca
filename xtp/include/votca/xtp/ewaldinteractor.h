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
#ifndef VOTCA_XTP_EWALDINTERACTOR_H
#define VOTCA_XTP_EWALDINTERACTOR_H

#include <array>
#include <boost/math/constants/constants.hpp>
#include <vector>

#include "bgsite.h"
#include "multipoleinteractor.h"

namespace votca {
namespace xtp {

class EwaldInteractor {
 public:
  EwaldInteractor(double alpha, double thole_damping)
      : mpInteract(alpha, thole_damping) {};

  // Fields
  Eigen::Vector3d r_staticFieldAtBy(const BGSite& site, const BGSite& nbSite,
                    const Eigen::Vector3d& shift);

 private:
  MultipoleInteractor mpInteract;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EANALYZE_H
