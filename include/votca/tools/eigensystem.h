/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_TOOLS_EIGENSYSTEM_H
#define __VOTCA_TOOLS_EIGENSYSTEM_H
#include "eigen.h"

namespace votca {
namespace tools {

class EigenSystem {
 public:
  // returns eigenvalues
  const Eigen::VectorXd& eigenvalues() const { return _eigenvalues; }
  Eigen::VectorXd& eigenvalues() { return _eigenvalues; }
  // returns eigenvectors
  const Eigen::MatrixXd& eigenvectors() const { return _eigenvectors; }
  Eigen::MatrixXd& eigenvectors() { return _eigenvectors; }
  // returns left eigenvectors or other in case of nonhermititan problem
  const Eigen::MatrixXd& eigenvectors2() const { return _eigenvectors_2; }
  Eigen::MatrixXd& eigenvectors2() { return _eigenvectors_2; }

  Eigen::ComputationInfo info() const { return _info; }
  Eigen::ComputationInfo& info() { return _info; }

  void clear() {
    _info = Eigen::Success;
    _eigenvalues.resize(0);
    _eigenvectors.resize(0, 0);
    _eigenvectors_2.resize(0, 0);
  }

 private:
  Eigen::ComputationInfo _info = Eigen::Success;
  Eigen::VectorXd _eigenvalues;
  Eigen::MatrixXd _eigenvectors;
  Eigen::MatrixXd _eigenvectors_2;
};

}  // namespace tools
}  // namespace votca
#endif  // __VOTCA_TOOLS_EIGENSYSTEM_H
