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
#ifndef VOTCA_XTP_ANDERSON_MIXING_H
#define VOTCA_XTP_ANDERSON_MIXING_H

// Standard includes
#include <vector>

// Local VOTCA includes
#include "votca/xtp/eigen.h"

namespace votca {
namespace xtp {
/**
 * \brief Anderson mixing as convergence acceleration in SCF/fixed point
 * problems
 *
 * Keeps a history of input and output solution vectors during a self-consistent
 * fixed-point procedure and determines a new solution from an optimized mixing.
 * Requires specification of the order of the method (maximum history to take)
 * and a mixing parameter.
 *
 *  B. Baumeier, Diploma Thesis, Appendix C (2005)
 *
 *  I. Ramiere and T. Helfer, Iterative Residual-Based Vector Methods to
 *  Accelerate Fixed Point Iterations, Comput. Math. Appl. 70, 2210 (2015)
 */
class Anderson {
 public:
  const Eigen::VectorXd MixHistory();

  void UpdateOutput(const Eigen::VectorXd &newOutput);
  void UpdateInput(const Eigen::VectorXd &newInput);
  void Configure(const Index order, const double alpha);

 private:
  std::vector<Eigen::VectorXd> _input;
  std::vector<Eigen::VectorXd> _output;

  double _alpha = 0.7;
  Index _order = 25;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ANDERSON_MIXING_H
