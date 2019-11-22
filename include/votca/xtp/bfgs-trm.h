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
#ifndef VOTCA_XTP_BFGSTRM_H
#define VOTCA_XTP_BFGSTRM_H

#include <functional>
#include <vector>
#include <votca/xtp/logger.h>
#include <votca/xtp/optimiser_costfunction.h>

namespace votca {
namespace xtp {

class BFGSTRM {
 public:
  BFGSTRM(Optimiser_costfunction& costfunction)
      : _costfunction(costfunction), _logging(false) {
    _hessian = Eigen::MatrixXd::Identity(costfunction.NumParameters(),
                                         costfunction.NumParameters());
  }

  void setLog(Logger* pLog) {
    _logging = true;
    _pLog = pLog;
  }

  void setTrustRadius(double trust_radius) { _trust_radius = trust_radius; }

  double getTrustRadius() const { return _trust_radius; }

  void setCallbacks(const std::vector<std::function<void()> >& callbacks) {
    _callbacks = callbacks;
  }

  void setNumofIterations(Index iterations) { _max_iteration = iterations; }

  void Optimize(const Eigen::VectorXd& initialparameters);

  bool Success() const { return _success; }
  std::string getErrorMessage() const { return _errormessage; }

  double getCost() const { return _cost; }

  Index getIteration() const { return _iteration; }

  const Eigen::VectorXd getParameters() const { return _parameters; }

  void setInitialHessian(const Eigen::MatrixXd& hessian) { _hessian = hessian; }

 private:
  Optimiser_costfunction& _costfunction;

  void UpdateHessian(const Eigen::VectorXd& delta_pos,
                     const Eigen::VectorXd& delta_gradient);
  double QuadraticEnergy(const Eigen::VectorXd& gradient,
                         const Eigen::VectorXd& delta_pos) const;
  bool AcceptRejectStep(const Eigen::VectorXd& delta_pos,
                        const Eigen::VectorXd& gradient, double energy_delta);

  std::string _errormessage;
  bool _success = true;
  bool _logging;
  Index _iteration = 0;

  std::vector<std::function<void()> > _callbacks;

  Eigen::MatrixXd _hessian;
  Eigen::VectorXd _parameters;

  double _cost = std::numeric_limits<double>::max();

  double _trust_radius = 0.1;

  Index _max_iteration = 200;

  Logger* _pLog;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_BFGSTRM_H
