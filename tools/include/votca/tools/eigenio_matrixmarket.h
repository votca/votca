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

#ifndef VOTCA_TOOLS_EIGENIO_MATRIXMARKET_H
#define VOTCA_TOOLS_EIGENIO_MATRIXMARKET_H

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace tools {

namespace EigenIO_MatrixMarket {

Eigen::VectorXd ReadVector(const std::string& filename);

void WriteVector(const std::string& filename, const Eigen::VectorXd& output);

void WriteMatrix(const std::string& filename, const Eigen::MatrixXd& output);

Eigen::MatrixXd ReadMatrix(const std::string& filename);

}  // namespace EigenIO_MatrixMarket

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_EIGENIO_MATRIXMARKET_H
