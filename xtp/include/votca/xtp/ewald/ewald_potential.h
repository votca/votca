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
#ifndef VOTCA_XTP_VXC_POTENTIAL_H
#define VOTCA_XTP_VXC_POTENTIAL_H

// Local VOTCA includes
#include "votca/xtp/grid_containers.h"
#include "votca/xtp/gridbox.h"

#undef LOG

namespace votca {
namespace xtp {

template <class Grid>
class Ewald_Potential {
 public:
  explicit Ewald_Potential(const Grid& grid) : grid_(grid) {};
  ~Ewald_Potential();

  Mat_p_Energy IntegrateEwald(const Eigen::MatrixXd& density_matrix) const;

 private:

  //XC_entry EvaluateXC(double rho, double sigma) const;

  const Grid grid_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_EWALD_POTENTIAL_H
