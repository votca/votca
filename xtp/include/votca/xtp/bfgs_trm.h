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
#ifndef VOTCA_XTP_BFGS_TRM_H
#define VOTCA_XTP_BFGS_TRM_H

// Standard includes
#include <functional>
#include <vector>

// Local VOTCA includes
#include "logger.h"
#include "optimiser_costfunction.h"

namespace votca {
namespace xtp {

class BFGSTRM {
 public:
  BFGSTRM(Optimiser_costfunction& costfunction) : costfunction_(costfunction) {
    hessian_ = Eigen::MatrixXd::Identity(costfunction.NumParameters(),
                                         costfunction.NumParameters());
  }

  void setLog(Logger* pLog) { pLog_ = pLog; }

  void setTrustRadius(double trust_radius) { trust_radius_ = trust_radius; }

  double getTrustRadius() const { return trust_radius_; }

  void setCallbacks(const std::vector<std::function<void()> >& callbacks) {
    callbacks_ = callbacks;
  }

  void setNumofIterations(Index iterations) { max_iteration_ = iterations; }

  void Optimize(const Eigen::VectorXd& initialparameters);

  bool Success() const { return success_; }
  std::string getErrorMessage() const { return errormessage_; }

  double getCost() const { return cost_; }

  Index getIteration() const { return iteration_; }

  const Eigen::VectorXd getParameters() const { return parameters_; }

  void setInitialHessian(const Eigen::MatrixXd& hessian) { hessian_ = hessian; }

 private:
  Optimiser_costfunction& costfunction_;

  void UpdateHessian(const Eigen::VectorXd& delta_pos,
                     const Eigen::VectorXd& delta_gradient);
  double QuadraticEnergy(const Eigen::VectorXd& gradient,
                         const Eigen::VectorXd& delta_pos) const;
  bool AcceptRejectStep(const Eigen::VectorXd& delta_pos,
                        const Eigen::VectorXd& gradient, double energy_delta);

  std::string errormessage_;
  bool success_ = true;
  Index iteration_ = 0;

  std::vector<std::function<void()> > callbacks_;

  Eigen::MatrixXd hessian_;
  Eigen::VectorXd parameters_;

  double cost_ = std::numeric_limits<double>::max();

  double trust_radius_ = 0.1;

  Index max_iteration_ = 200;

  Logger* pLog_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_BFGS_TRM_H
