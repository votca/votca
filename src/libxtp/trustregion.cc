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

#include <votca/xtp/trustregion.h>
#include <iostream>
#include <vector>
#include <iomanip>

namespace votca {
  namespace xtp {

    Eigen::VectorXd TrustRegion::CalculateStep(const Eigen::VectorXd& gradient, const Eigen::MatrixXd& Hessian, double delta)const {
      //calculate unrestricted step
      Eigen::VectorXd freestep = Hessian.colPivHouseholderQr().solve(-gradient);

      //if inside use the step;
      if (freestep.norm() < delta) {
        return freestep;
      }

      //calculate step on the boundary
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hessian);
      const Eigen::VectorXd factor = (es.eigenvectors().transpose() * gradient).cwiseAbs2();
      double lambda = 0;
      //hard case
      if (std::abs(factor[0]) < 1e-18) {
        std::vector<int> index;
        Eigen::VectorXd diag = Eigen::VectorXd::Zero(factor.size());
        //construct pseudo inverse
        for (int i = 0; i < factor.size(); i++) {
          double entry = es.eigenvalues()(i) - es.eigenvalues()(0);
          if (entry < 1e-14) {
            index.push_back(i);
            continue;
          } else {
            diag(i) = 1 / entry;
          }
        }
        Eigen::MatrixXd pseudo_inv = es.eigenvectors() * diag.asDiagonal() * es.eigenvectors().transpose();
        Eigen::VectorXd step = -pseudo_inv*gradient;
        if (step.norm() <= delta) {
          double tau = std::sqrt(delta * delta - step.squaredNorm());
          step += tau * es.eigenvectors().col(index[0]);
          return step;
        }
      }

      //sort out all the factors for the small eigenvalues, 
      int start_index = 0;
      for (; start_index < factor.size(); start_index++) {
        if (factor[start_index] < 1e-18) {
          continue;
        }
        break;
      }

      // start value for lambda  a bit higher than lowest eigenvalue of Hessian
      lambda = -es.eigenvalues()(start_index) + std::sqrt(factor(start_index)) / delta;
      TrustRegionFunction TRF = TrustRegionFunction(factor, es, delta, start_index);
      for (int iter = 0; iter < 100; iter++) {
        std::pair<double, double> result = TRF.Evaluate(lambda);
        double func_value = result.first;
        double update = result.second;
        lambda += update;
        if (update < 1e-14 || std::abs(func_value) < 1e-12) {
          break;
        }
      }

      //this is effectively the solution of (H+I\lambda)*\Delta p=-g with \lambda adjusted so that ||p||=delta 
      Eigen::VectorXd new_delta_pos = Eigen::VectorXd::Zero(gradient.size());
      for (int i = start_index; i < gradient.size(); i++) {
        new_delta_pos -= es.eigenvectors().col(i) * (es.eigenvectors().col(i).transpose() * gradient) / (es.eigenvalues()(i) + lambda);
      }
      return new_delta_pos;
    }

  }
}
