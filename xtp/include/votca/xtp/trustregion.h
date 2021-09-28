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
#ifndef VOTCA_XTP_TRUSTREGION_H
#define VOTCA_XTP_TRUSTREGION_H

// Standard includes
#include <iostream>

// Local VOTCA includes
#include "eigen.h"

// Solves the trustregion subproblem g^T*s+0.5*s^T H s = min with ||s||<=delta

namespace votca {
namespace xtp {

class TrustRegion {
 public:
  Eigen::VectorXd CalculateStep(const Eigen::VectorXd& gradient,
                                const Eigen::MatrixXd& Hessian,
                                double delta) const;

 private:
  class TrustRegionFunction {
   public:
    TrustRegionFunction(
        const Eigen::VectorXd& factor,
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& hessian,
        double trust_radius, Index startindex)
        : factor_(factor),
          hessian_(hessian),
          trust_radius_(trust_radius),
          startindex_(startindex) {
      ;
    }

    // Calculates \phi and \phi/\phi'
    std::pair<double, double> Evaluate(double lambda) {
      Index size = factor_.size() - startindex_;
      Eigen::ArrayXd quotient =
          (hessian_.eigenvalues().array() + lambda).tail(size);
      const double p2 = (factor_.array().tail(size) / (quotient.pow(2))).sum();
      const double p = std::sqrt(p2);
      const double q2 = (factor_.array().tail(size) / (quotient.pow(3))).sum();
      std::pair<double, double> result;
      result.first = 1 / trust_radius_ - 1 / p;
      result.second = p2 / q2 * (p - trust_radius_) / trust_radius_;
      return result;
    }

   private:
    const Eigen::VectorXd& factor_;
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& hessian_;
    double trust_radius_;
    Index startindex_;
  };
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_TRUSTREGION_H
