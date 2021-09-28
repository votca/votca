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

#ifndef VOTCA_XTP_GAUSS_LAGUERRE_QUADRATURE_H
#define VOTCA_XTP_GAUSS_LAGUERRE_QUADRATURE_H

#include "votca/xtp/GaussianQuadratureBase.h"

namespace votca {
namespace xtp {

class Gauss_Laguerre_Quadrature : public GaussianQuadratureBase {

 public:
  double ScaledPoint(Index i) const final { return points_[i]; }

  double ScaledWeight(Index i) const final { return weights_[i]; }

 protected:
  // The laguerre method is suitable for integration limits a = 0 b =
  // +infty. Here we have value1 and value2 because we split the original
  // integral from -infty to +infty in two parts
  bool UseSymmetry() const final { return true; }
  void FillPoints() final;
  void FillAdaptedWeights() final;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GAUSS_LAGUERRE_QUADRATURE_H
