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
#ifndef VOTCA_XTP_ADIIS_COSTFUNCTION_H
#define VOTCA_XTP_ADIIS_COSTFUNCTION_H

// Local VOTCA includes
#include "optimiser_costfunction.h"

namespace votca {
namespace xtp {

class ADIIS_costfunction : public Optimiser_costfunction {
 public:
  ADIIS_costfunction(Eigen::VectorXd DiF, Eigen::MatrixXd DiFj) {
    DiF_ = DiF;
    DiFj_ = DiFj;
  }

  double EvaluateCost(const Eigen::VectorXd& parameters) override {
    Eigen::VectorXd c = parameters.cwiseAbs2();
    double xnorm = c.sum();
    c /= xnorm;
    return (2 * c.transpose() * DiF_ + c.transpose() * DiFj_ * c).value();
  }

  Eigen::VectorXd EvaluateGradient(const Eigen::VectorXd& parameters) override {
    Eigen::VectorXd c = parameters.cwiseAbs2();
    double xnorm = c.sum();
    c /= xnorm;
    Eigen::VectorXd dEdc = 2.0 * DiF_ + DiFj_ * c + DiFj_.transpose() * c;
    Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(c.size(), c.size());
    for (Index i = 0; i < jac.rows(); i++) {
      for (Index j = 0; j < jac.cols(); j++) {
        jac(i, j) = -c(i) * 2.0 * parameters(j) / xnorm;
      }
      // Extra term on diagonal
      jac(i, i) += 2.0 * parameters(i) / xnorm;
    }
    return jac.transpose() * dEdc;
  }

  Index NumParameters() const override { return Index(DiF_.size()); }

  bool Converged(const Eigen::VectorXd&, double,
                 const Eigen::VectorXd& gradient) override {
    return gradient.cwiseAbs().maxCoeff() < 1.e-7;
  }

 private:
  Eigen::VectorXd DiF_;
  Eigen::MatrixXd DiFj_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_ADIIS_COSTFUNCTION_H
