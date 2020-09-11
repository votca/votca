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

#ifndef VOTCA_XTP_GAUSS_LEGENDRE_QUADRATURE_H
#define VOTCA_XTP_GAUSS_LEGENDRE_QUADRATURE_H

#include "votca/xtp/GaussianQuadratureBase.h"
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

class Gauss_Legendre_Quadrature_Base : public GaussianQuadratureBase {

 protected:
  void FillPoints() final;
  void FillAdaptedWeights() final;
};

class Gauss_modified_Legendre_Quadrature
    : public Gauss_Legendre_Quadrature_Base {
 public:
  double ScaledPoint(Index i) const final {
    double exponent = (1.0 + points_(i)) / (1.0 - points_(i));
    return std::pow(0.5, exponent);
  }
  double ScaledWeight(Index i) const final {
    double den = (1.0 - weights_(i)) * (1.0 - weights_(i));
    return (2.0 * weights_(i) * 0.5) / den;
  }

 protected:
  // The modified legendre method is suitable for integration limits a = 0 b =
  // +infty. Here we have value1 and value2 because we split the original
  // integral from -infty to +infty in two parts. Original legendre quadrature
  // is meant for integratal with integration limits of -1 and 1. To overcome
  // this we use the transformation x' = 0.5 ^ (1+x/1-x)
  bool UseSymmetry() const final { return true; }
};

class Gauss_Legendre_Quadrature : public Gauss_Legendre_Quadrature_Base {
 public:
  double ScaledPoint(Index i) const final {
    return std::tan(0.5 * votca::tools::conv::Pi * points_(i));
  }

  double ScaledWeight(Index i) const final {
    const double halfpi = 0.5 * tools::conv::Pi;
    double den = std::cos(halfpi * points_[i]) * std::cos(halfpi * points_[i]);
    return weights_(i) * halfpi / den;
  }

 protected:
  // This particular legendre quadrature is suitable for integration limits a
  // = -infty b = +infty. Original Legendre quadrature is meant for
  // integration limits -1 and +1. The change of variables is x' = tan (pi/2 *
  // x). When x=-1 we have x'=-infty. When x=1 we have x'=+infty
  bool UseSymmetry() const final { return false; }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GAUSS_LEGENDRE_QUADRATURE_H