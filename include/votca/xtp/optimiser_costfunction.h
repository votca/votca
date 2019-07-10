/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#ifndef VOTCA_XTP_OPTIMISER_COSTFUNCTION_H
#define VOTCA_XTP_OPTIMISER_COSTFUNCTION_H

#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

class Optimiser_costfunction {
 public:
  virtual ~Optimiser_costfunction(){};

  virtual double EvaluateCost(const Eigen::VectorXd& parameters) = 0;

  virtual Eigen::VectorXd EvaluateGradient(
      const Eigen::VectorXd& parameters) = 0;

  virtual int NumParameters() const = 0;

  virtual bool Converged(const Eigen::VectorXd& delta_parameters,
                         double delta_cost,
                         const Eigen::VectorXd& gradient) = 0;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_OPTIMISER_COSTFUNCTION_H
