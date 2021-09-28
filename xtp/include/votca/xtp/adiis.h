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
#ifndef VOTCA_XTP_ADIIS_H
#define VOTCA_XTP_ADIIS_H

// Standard includes
#include <memory>
#include <vector>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

class ADIIS {
 public:
  Eigen::VectorXd CalcCoeff(const std::vector<Eigen::MatrixXd>& dmathist,
                            const std::vector<Eigen::MatrixXd>& mathist);

  bool Info() { return success; }

 private:
  bool success = true;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ADIIS_H
