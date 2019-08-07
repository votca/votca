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
#ifndef __XTP_TRUSTREGION__H
#define __XTP_TRUSTREGION__H

#include <iostream>
#include <votca/xtp/eigen.h>

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
        double trust_radius, int startindex)
        : _factor(factor),
          _hessian(hessian),
          _trust_radius(trust_radius),
          _startindex(startindex) {
      ;
    }

    // Calculates \phi and \phi/\phi'
    std::pair<double, double> Evaluate(double lambda) {
      int size = _factor.size() - _startindex;
      Eigen::ArrayXd quotient =
          (_hessian.eigenvalues().array() + lambda).tail(size);
      const double p2 = (_factor.array().tail(size) / (quotient.pow(2))).sum();
      const double p = std::sqrt(p2);
      const double q2 = (_factor.array().tail(size) / (quotient.pow(3))).sum();
      std::pair<double, double> result;
      result.first = 1 / _trust_radius - 1 / p;
      result.second = p2 / q2 * (p - _trust_radius) / _trust_radius;
      return result;
    }

   private:
    const Eigen::VectorXd& _factor;
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& _hessian;
    double _trust_radius;
    int _startindex;
  };
};

}  // namespace xtp
}  // namespace votca
#endif /* TRUSTREGION_H */
