
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
#ifndef VOTCA_XTP_MPINTERACTIONS_H
#define VOTCA_XTP_MPINTERACTIONS_H

#include <algorithm>
#include <vector>
// Local VOTCA includes
#include "votca/xtp/interactiontensor.h"
#include "votca/xtp/multipole.h"

namespace votca {
namespace xtp {

class MPInteractions {
 public:
 MPInteractions(double alpha, double thole_damping);

 Eigen::Vector3d r_fieldAtBy(const Multipole& mp2, const Eigen::Vector3d& dr);
  

 private:
 double alpha_;
 double thole;
 double thole2;
 double thole3;

};  

}  // namespace xtp
}  // namespace votca

#endif