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

/**
 * Pulay direct inversion in the iterative subspace.
 *
 * The class stores recent commutator error matrices and solves the augmented
 * linear DIIS system to obtain extrapolation coefficients for the next Fock
 * guess.
 */
class DIIS {
 public:
  /// Store a new restricted-spin DIIS error matrix and update the B matrix
  /// history.
  void Update(Index maxerrorindex, const Eigen::MatrixXd& errormatrix);
  /// Store new alpha and beta DIIS error matrices and update the
  /// unrestricted-spin B matrix history.
  void Update(Index maxerrorindex, const Eigen::MatrixXd& errormatrix_alpha,
              const Eigen::MatrixXd& errormatrix_beta);
  /// Solve the Pulay linear system and return the DIIS interpolation
  /// coefficients.
  Eigen::VectorXd CalcCoeff();

  /// Set the maximum number of past error matrices retained for DIIS
  /// extrapolation.
  void setHistLength(Index length) { histlength_ = length; }

  /// Report whether the most recent DIIS setup and solve completed
  /// successfully.
  bool Info() { return success; }

 private:
  bool success = true;
  Index histlength_;
  std::vector<std::vector<double> > Diis_Bs_;
  std::vector<Eigen::MatrixXd> errormatrixhist_;
  std::vector<Eigen::MatrixXd> errormatrixhist_alpha_;
  std::vector<Eigen::MatrixXd> errormatrixhist_beta_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DIIS_H
