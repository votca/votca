/*
 *            Copyright 2009-2022 The VOTCA Development Team
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
#ifndef VOTCA_XTP_PMLOCALIZATION_H
#define VOTCA_XTP_PMLOCALIZATION_H
#include "logger.h"
#include "votca/tools/property.h"
#include "votca/xtp/orbitals.h"

namespace votca {
namespace xtp {

class PMLocalization {
 public:
  PMLocalization(Logger &log, const tools::Property &options) : log_(log) {
    nrOfIterations_ = options.get(".max_iterations").as<Index>();
    convergence_limit_ = options.get(".convergence_limit").as<double>();
  };
  void computePML(Orbitals &orbitals);

 private:
  Logger &log_;
  Eigen::MatrixX2d rotateorbitals(const Eigen::MatrixX2d &maxorbs, Index s,
                                  Index t);

  void initial_penalty();
  void update_penalty(Index s, Index t);
  void check_orthonormality();
  Eigen::VectorXd calculate_lmo_energies(Orbitals &orbitals);
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> sort_lmos(
      const Eigen::VectorXd &energies);

  Eigen::VectorXd pop_per_atom(const Eigen::VectorXd &orbital);
  Eigen::Vector2d offdiag_penalty_elements(const Eigen::MatrixXd &s_overlap,
                                           Index s, Index t);

  Eigen::MatrixXd localized_orbitals_;

  AOBasis aobasis_;
  Eigen::MatrixXd overlap_;
  Eigen::MatrixXd A_;
  Eigen::MatrixXd B_;

  Eigen::MatrixXd PM_penalty_;

  Eigen::MatrixXd MullikenPop_orb_per_atom_;

  std::vector<Index> numfuncpatom_;

  Index nrOfIterations_ = 0;
  double convergence_limit_ = 0.0;
};

}  // namespace xtp
}  // namespace votca
#endif