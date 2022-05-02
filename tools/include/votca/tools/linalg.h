/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_TOOLS_LINALG_H
#define VOTCA_TOOLS_LINALG_H

// Local VOTCA includes
#include "eigen.h"
#include "eigensystem.h"
#include "types.h"

namespace votca {
namespace tools {

/**
 * \brief solves A*x=b under the constraint B*x = 0
 * @return x
 * @param A matrix for linear equation system
 * @param b inhomogenity
 * @param constr constrained condition
 *
 * This function implements the qrsolver under constraints
 */
Eigen::VectorXd linalg_constrained_qrsolve(const Eigen::MatrixXd& A,
                                           const Eigen::VectorXd& b,
                                           const Eigen::MatrixXd& constr);

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_LINALG_H
