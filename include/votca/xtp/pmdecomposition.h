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
#ifndef VOTCA_XTP_PMDECOMPOSITION_H
#define VOTCA_XTP_PMDECOMPOSITION_H
#include "logger.h"
#include "votca/xtp/orbitals.h"

namespace votca {
namespace xtp {

class PMDecomposition {
 public:
  PMDecomposition() = default;
  PMDecomposition(Orbitals &orbitals_, Logger &log_)
      : orbitals(orbitals_), log(log_){};
  void compute();

 private:
  Orbitals orbitals;
  Logger &log;
  Eigen::MatrixXd rotatedorbitals(Eigen::MatrixXd &maxorbs,
                                  Index s, Index t);
  Eigen::MatrixXd orbitalselections(Eigen::MatrixXd &m,
                                    const Eigen::MatrixXd &S);
  void update_maximums(Eigen::MatrixXd &m, Index col1, Index col2,
                       Eigen::MatrixXd &new_orbs);
Eigen::MatrixXd columnwise(const Eigen::MatrixXd &S, Eigen::VectorXd &v);
Eigen::MatrixXd rowwise(const Eigen::MatrixXd &S, Eigen::VectorXd &v);
BasisSet basis;
AOBasis aobasis;
Eigen::MatrixXd A;
Eigen::MatrixXd B;
};
}  // namespace xtp
}  // namespace votca
#endif