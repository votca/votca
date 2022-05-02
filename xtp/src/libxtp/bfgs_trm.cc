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

// Third party includes
#include <boost/format.hpp>

// Local VOTCA includes
#include "votca/xtp/atom.h"
#include "votca/xtp/bfgs_trm.h"
#include "votca/xtp/trustregion.h"

namespace votca {
namespace xtp {

void BFGSTRM::Optimize(const Eigen::VectorXd& initialparameters) {
  parameters_ = initialparameters;
  cost_ = costfunction_.EvaluateCost(parameters_);
  double lastcost = cost_;
  Eigen::VectorXd gradient = costfunction_.EvaluateGradient(parameters_);
  for (auto& func : callbacks_) {
    func();
  }

  Eigen::VectorXd delta_p_trial = Eigen::VectorXd::Zero(parameters_.size());
  Eigen::VectorXd last_gradient = Eigen::VectorXd::Zero(parameters_.size());
  double delta_cost = 0;
  for (iteration_ = 1; iteration_ <= max_iteration_; iteration_++) {
    for (Index i = 0; i < 100; i++) {
      TrustRegion subproblem;
      delta_p_trial =
          subproblem.CalculateStep(gradient, hessian_, trust_radius_);
      double trialcost =
          costfunction_.EvaluateCost(parameters_ + delta_p_trial);
      delta_cost = trialcost - lastcost;
      bool step_accepted =
          AcceptRejectStep(delta_p_trial, gradient, delta_cost);
      if (step_accepted) {
        cost_ = trialcost;
        parameters_ += delta_p_trial;
        break;
      }
    }
    gradient = costfunction_.EvaluateGradient(parameters_);
    if (iteration_ > 1) {
      UpdateHessian(delta_p_trial, gradient - last_gradient);
    }
    lastcost = cost_;
    last_gradient = gradient;
    for (auto& func : callbacks_) {
      func();
    }
    if (costfunction_.Converged(delta_p_trial, delta_cost, gradient)) {
      break;
    } else if (iteration_ == max_iteration_) {
      success_ = false;
      XTP_LOG(Log::warning, *pLog_)
          << (boost::format("BFGS-TRM @iteration %1$d: not converged after "
                            "%2$d iterations ") %
              iteration_ % max_iteration_)
                 .str()
          << std::flush;
    }
  }
  return;
}

/* Accept/reject the new geometry and adjust trust radius, if required */
bool BFGSTRM::AcceptRejectStep(const Eigen::VectorXd& delta_p,
                               const Eigen::VectorXd& gradient,
                               double cost_delta) {
  bool step_accepted = false;
  if (cost_delta > 0.0) {
    // total energy has unexpectedly increased, half the trust radius
    trust_radius_ = 0.25 * trust_radius_;
    XTP_LOG(Log::warning, *pLog_)
        << (boost::format("BFGS-TRM @iteration %1$d: DeltaCost %2$2.4e step "
                          "rejected ") %
            iteration_ % cost_delta)
               .str()
        << std::flush;
    XTP_LOG(Log::warning, *pLog_)
        << (boost::format(
                "BFGS-TRM @iteration %1$d: new trust radius %2$2.4e") %
            iteration_ % trust_radius_)
               .str()
        << std::flush;

  } else {
    // total energy has decreased, we accept the step but might update the trust
    // radius
    step_accepted = true;
    // adjust trust radius, if required
    double tr_check = cost_delta / QuadraticEnergy(gradient, delta_p);
    double norm_delta_p = delta_p.squaredNorm();
    if (tr_check > 0.75 &&
        1.25 * norm_delta_p > trust_radius_ * trust_radius_) {
      trust_radius_ = 2.0 * trust_radius_;
    } else if (tr_check < 0.25) {
      trust_radius_ = 0.25 * trust_radius_;
    }
    XTP_LOG(Log::warning, *pLog_)
        << (boost::format(
                "BFGS-TRM @iteration %1$d: DeltaCost/QuadraticApprox %2$2.4f "
                "step accepted ") %
            iteration_ % tr_check)
               .str()
        << std::flush;
    XTP_LOG(Log::warning, *pLog_)
        << (boost::format(
                "BFGS-TRM @iteration %1$d: new trust radius %2$2.4e") %
            iteration_ % trust_radius_)
               .str()
        << std::flush;
  }
  return step_accepted;
}

void BFGSTRM::UpdateHessian(const Eigen::VectorXd& delta_pos,
                            const Eigen::VectorXd& delta_gradient) {
  // second term in BFGS update (needs current Hessian)
  hessian_ -= hessian_ * delta_pos * delta_pos.transpose() *
              hessian_.transpose() /
              (delta_pos.transpose() * hessian_ * delta_pos).value();
  // first term in BFGS update
  hessian_ += (delta_gradient * delta_gradient.transpose()) /
              (delta_gradient.transpose() * delta_pos);
  // symmetrize Hessian (since d2E/dxidxj should be symmetric)
  hessian_ = 0.5 * (hessian_ + hessian_.transpose());
  return;
}

/* Estimate energy change based on quadratic approximation */
double BFGSTRM::QuadraticEnergy(const Eigen::VectorXd& gradient,
                                const Eigen::VectorXd& delta_pos) const {
  return (gradient.transpose() * delta_pos).value() +
         0.5 * (delta_pos.transpose() * hessian_ * delta_pos).value();
}

}  // namespace xtp
}  // namespace votca
