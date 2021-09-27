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
#ifndef VOTCA_XTP_DIIS_H
#define VOTCA_XTP_DIIS_H

// Standard includes
#include <vector>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

class DIIS {
 public:
  void Update(Index maxerrorindex, const Eigen::MatrixXd& errormatrix);
  Eigen::VectorXd CalcCoeff();

  void setHistLength(Index length) { histlength_ = length; }

  bool Info() { return success; }

 private:
  bool success = true;
  Index histlength_;
  std::vector<std::vector<double> > Diis_Bs_;
  std::vector<Eigen::MatrixXd> errormatrixhist_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DIIS_H
