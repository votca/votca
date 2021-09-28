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
#ifndef VOTCA_XTP_NEWTON_RAPSON_H
#define VOTCA_XTP_NEWTON_RAPSON_H

// VOTCA includes
#include <votca/tools/types.h>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Newton Rapson rootfinder for 1d functions
 *
 * https://en.wikipedia.org/wiki/Newton%27s_method#Modified_Newton_methods
 *
 */
template <class Func>
class NewtonRapson {
 public:
  enum Errors { success, smalldenom, notconverged };
  NewtonRapson(Index max_iterations, double tolerance)
      : max_iterations_(max_iterations), tolerance_(tolerance) {}

  NewtonRapson(Index max_iterations, double tolerance, double alpha)
      : max_iterations_(max_iterations), tolerance_(tolerance), alpha_(alpha) {}

  double FindRoot(const Func& f, double x0) {
    info_ = Errors::notconverged;
    double x = x0;
    for (iter_ = 0; iter_ < max_iterations_; iter_++) {

      std::pair<double, double> res = f(x);
      if (std::abs(res.second) < 1e-12) {
        info_ = Errors::smalldenom;
        break;
      }

      double step = -alpha_ * res.first / res.second;
      if (std::abs(step) < tolerance_) {
        info_ = Errors::success;
        break;
      }

      x += step;
    }

    return x;
  }

  Errors getInfo() const { return info_; }
  Index getIterations() const { return iter_; }

 private:
  Errors info_ = Errors::notconverged;
  Index max_iterations_;
  Index iter_;
  double tolerance_;
  double alpha_ = 1.0;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_NEWTON_RAPSON_H
