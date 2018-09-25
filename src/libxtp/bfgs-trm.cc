/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <votca/xtp/bfgs-trm.h>
#include <boost/format.hpp>
#include <votca/xtp/trustregion.h>

namespace votca {
  namespace xtp {

    void BFGSTRM::Optimize(const Eigen::VectorXd& initialparameters) {
      _parameters = initialparameters;

      Eigen::VectorXd gradient = Eigen::VectorXd::Zero(_parameters.size());
      Eigen::VectorXd delta_p_trial = Eigen::VectorXd::Zero(_parameters.size());

      double lastcost = _costfunction.EvaluateCost(_parameters);
      Eigen::VectorXd last_gradient = Eigen::VectorXd::Zero(_parameters.size());
      double delta_cost = 0;
      std::cout << "here1" << std::endl;
      for (_iteration = 0; _iteration < _max_iteration; _iteration++) {
        gradient = _costfunction.EvaluateGradient(_parameters);
        std::cout << "here3" << std::endl;
        bool step_accepted = false;
        for (int i = 0; i < 100; i++) {
          TrustRegion subproblem;
          delta_p_trial = subproblem.CalculateStep(gradient, _hessian, _trust_radius);
          std::cout << "here2" << std::endl;
          double trialcost = _costfunction.EvaluateCost(_parameters + delta_p_trial);
          delta_cost = trialcost - lastcost;
          step_accepted = AcceptRejectStep(delta_p_trial, gradient, delta_cost);
          std::cout << "here4" << std::endl;
          if (step_accepted) {
            _cost = trialcost;
            _parameters += delta_p_trial;
            break;
          }
          std::cout << "here4a" << std::endl;
        }
        std::cout << "here5" << std::endl;
        if (_iteration > 0) {
          UpdateHessian(delta_p_trial, gradient - last_gradient);
        }
        std::cout << "here6" << std::endl;
        lastcost = _cost;
        last_gradient = gradient;
        for (auto& func : _callbacks) {
          func();
        }
        if (_costfunction.Converged(delta_p_trial, delta_cost, gradient)) {
          break;
        } else if (_iteration == _max_iteration - 1) {
          _success = false;
          if (_logging) {
            CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: not converged after %2$d iterations ")
                    % _iteration % _max_iteration).str() << std::flush;
          }
        }
        std::cout << "here6" << std::endl;
      }
      return;
    }

    /* Accept/reject the new geometry and adjust trust radius, if required */
    bool BFGSTRM::AcceptRejectStep(const Eigen::VectorXd& delta_p, const Eigen::VectorXd& gradient, double cost_delta) {
      bool step_accepted = false;
      if (cost_delta > 0.0) {
        // total energy has unexpectedly increased, half the trust radius
        _trust_radius = 0.25 * _trust_radius;
        if (_logging) {
          CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: step rejected ")
                  % _iteration).str() << std::flush;
          CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: new trust radius %2$8.10e")
                  % _iteration % _trust_radius).str() << std::flush;
        }
      } else {
        // total energy has decreased, we accept the step but might update the trust radius
        step_accepted = true;
        // adjust trust radius, if required
        double tr_check = cost_delta / QuadraticEnergy(gradient, delta_p);
        double norm_delta_p = delta_p.squaredNorm();
        if (tr_check > 0.75 && 1.25 * norm_delta_p > _trust_radius * _trust_radius) {
          _trust_radius = 2.0 * _trust_radius;
        } else if (tr_check < 0.25) {
          _trust_radius = 0.25 * _trust_radius;
        }
        if (_logging) {
          CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: step accepted ")
                  % _iteration).str() << std::flush;
          CTP_LOG(ctp::logINFO, *_pLog) << (boost::format("BFGS-TRM @iteration %1$d: new trust radius %2$8.10e")
                  % _iteration % _trust_radius).str() << std::flush;
        }
      }
      return step_accepted;
    }

    void BFGSTRM::UpdateHessian(const Eigen::VectorXd& delta_pos, const Eigen::VectorXd& delta_gradient) {
      // second term in BFGS update (needs current Hessian)
      _hessian -= _hessian * delta_pos * delta_pos.transpose() * _hessian.transpose() / (delta_pos.transpose() * _hessian * delta_pos).value();
      // first term in BFGS update
      _hessian += (delta_gradient * delta_gradient.transpose()) / (delta_gradient.transpose() * delta_pos);
      // symmetrize Hessian (since d2E/dxidxj should be symmetric)
      _hessian = 0.5 * (_hessian + _hessian.transpose());
      return;
    }

    /* Estimate energy change based on quadratic approximation */
    double BFGSTRM::QuadraticEnergy(const Eigen::VectorXd& gradient, const Eigen::VectorXd& delta_pos) const {
      return (gradient.transpose() * delta_pos).value() + 0.5 * (delta_pos.transpose() * _hessian * delta_pos).value();
    }

  }
}
