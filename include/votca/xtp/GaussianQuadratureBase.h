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

#ifndef VOTCA_XTP_GAUSSIANQUADRATUREBASE_H
#define VOTCA_XTP_GAUSSIANQUADRATUREBASE_H

// Local VOTCA includes
#include "eigen.h"
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

class GaussianQuadratureBase {
 public:
  GaussianQuadratureBase() = default;

  void configure(Index order) {
    FillPoints();
    FillAdaptedWeights();
    CheckOrder(order, _map_points);
    CheckOrder(order, _map_AdaptedWeights);
    points_ = _map_points[order];
    weights_ = _map_AdaptedWeights[order];
  }

  virtual ~GaussianQuadratureBase() = default;

  Index Order() const { return points_.size(); }

  template <typename F>
  double Integrate(const F& f) const {
    double result = 0.0;
    for (Index j = 0; j < Order(); ++j) {
      result += ScaledWeight(j) * f(j, ScaledPoint(j), UseSymmetry());
    }
    return 0.5 * result / tools::conv::Pi;
  }
  virtual double ScaledPoint(Index i) const = 0;

  virtual double ScaledWeight(Index i) const = 0;

 private:
  void CheckOrder(Index order,
                  const std::map<Index, Eigen::VectorXd>& map) const;

 protected:
  Eigen::VectorXd points_;
  Eigen::VectorXd weights_;

  std::map<Index, Eigen::VectorXd> _map_points;
  std::map<Index, Eigen::VectorXd> _map_AdaptedWeights;

  virtual bool UseSymmetry() const = 0;

  virtual void FillPoints() = 0;
  virtual void FillAdaptedWeights() = 0;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GAUSSIANQUADRATUREBASE_H
