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

#ifndef VOTCA_XTP_SIGMA_EXACT_H
#define VOTCA_XTP_SIGMA_EXACT_H

// Local VOTCA includes
#include "rpa.h"
#include "sigma_base.h"

namespace votca {
namespace xtp {

class TCMatrix_gwbse;
class RPA;

class Sigma_Exact : public Sigma_base {

 public:
  Sigma_Exact(TCMatrix_gwbse& Mmn, RPA& rpa) : Sigma_base(Mmn, rpa){};

  // Sets up the screening parametrisation
  void PrepareScreening() final;
  // Calculates Sigma_c diagonal elements
  double CalcCorrelationDiagElement(Index gw_level,
                                    double frequency) const final;

  double CalcCorrelationDiagElementDerivative(Index gw_level,
                                              double frequency) const final;
  // Calculates Sigma_c off-diagonal elements
  double CalcCorrelationOffDiagElement(Index gw_level1, Index gw_level2,
                                       double frequency1,
                                       double frequency2) const final;

 private:
  RPA::rpa_eigensolution _rpa_solution;    // Eigenvalues, eigenvectors from RPA
  std::vector<Eigen::MatrixXd> _residues;  // Residues

  Eigen::MatrixXd CalcResidues(Index gw_level) const;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SIGMA_EXACT_H
