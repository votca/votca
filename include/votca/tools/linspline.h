/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _LINSPLINE_H
#define _LINSPLINE_H

#include "spline.h"
#include <iostream>
#include <votca/tools/eigen.h>

namespace votca {
namespace tools {

/**
 * \brief A Linear Spline Class
 *
 *  class supports linear interpolation and linear fit of data
 */

class LinSpline : public Spline {
 public:
  // default constructor
  LinSpline() = default;
  // LinSpline() :
  //    _boundaries(splineNormal) {}

  // destructor
  ~LinSpline() override = default;

  // construct an interpolation spline
  // x, y are the the points to construct interpolation, both vectors must be of
  // same size
  void Interpolate(const Eigen::VectorXd &x, const Eigen::VectorXd &y) override;

  // fit spline through noisy data
  // x,y are arrays with noisy data, both vectors must be of same size
  void Fit(const Eigen::VectorXd &x, const Eigen::VectorXd &y) override;

  // Calculate the function value
  double Calculate(double x) override;

  // Calculate the function derivative
  double CalculateDerivative(double x) override;
  using Spline::Calculate;
  using Spline::CalculateDerivative;

 protected:
  // a,b for piecewise splines: ax+b
  Eigen::VectorXd a;
  Eigen::VectorXd b;
};

inline double LinSpline::Calculate(double r) {
  Index interval = getInterval(r);
  return a(interval) * r + b(interval);
}

inline double LinSpline::CalculateDerivative(double r) {
  Index interval = getInterval(r);
  return a(interval);
}

}  // namespace tools
}  // namespace votca

#endif /* _LINSPLINE_H */
