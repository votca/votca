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

#include <iostream>
#include <votca/tools/akimaspline.h>
#include <votca/tools/linalg.h>

namespace votca {
namespace tools {

using namespace std;

void AkimaSpline::Interpolate(Eigen::VectorXd &x, Eigen::VectorXd &y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument(
        "error in AkimaSpline::Interpolate : sizes of vectors x and y do not "
        "match");
  }

  // Akima splines require at least 4 points
  if (x.size() < 4) {
    throw std::invalid_argument(
        "error in AkimaSpline::Interpolate : vectors x and y have to contain "
        "at least 4 points");
  }

  const Index N = x.size();

  // copy the grid points into f
  _r = x;

  // initialize vectors p1,p2,p3,p4 and t
  p0 = Eigen::VectorXd::Zero(N);
  p1 = Eigen::VectorXd::Zero(N);
  p2 = Eigen::VectorXd::Zero(N);
  p3 = Eigen::VectorXd::Zero(N);
  t = Eigen::VectorXd::Zero(N);

  double m1, m2, m3, m4;

  double temp, g0, g1, g2, x1, x2, y1, y2, x4, x5, y4, y5, m5;

  // boundary conditions
  // >> determine t(0), t(1) and t(N-2), t(N-1)
  switch (_boundaries) {
    case splineNormal:
      // Akima method: estimation of two more points on each side using a
      // degree two polyomial
      // resulting slopes t(0), t(1), t(N-2), t(N-1) are directly calculated

      // left side: t(0), t(1)
      temp = (x(1) - x(0)) / (x(2) - x(0));
      temp = temp * temp;
      g0 = y(0);
      g1 = ((y(1) - y(0)) - temp * (y(2) - y(0))) /
           ((x(1) - x(0)) - temp * (x(2) - x(0)));
      g2 = ((y(2) - y(0)) - g1 * (x(2) - x(0))) /
           ((x(2) - x(0)) * (x(2) - x(0)));
      x1 = x(0) - (x(2) - x(0));
      x2 = x(1) - (x(2) - x(0));
      y1 = g0 + g1 * (x1 - x(0)) + g2 * (x1 - x(0)) * (x1 - x(0));
      y2 = g0 + g1 * (x2 - x(0)) + g2 * (x2 - x(0)) * (x2 - x(0));
      m1 = (y2 - y1) / (x2 - x1);
      m2 = (y(0) - y2) / (x(0) - x2);
      m3 = (y(1) - y(0)) / (x(1) - x(0));
      m4 = (y(2) - y(1)) / (x(2) - x(1));
      t(0) = getSlope(m1, m2, m3, m4);
      m5 = (y(3) - y(2)) / (x(3) - x(2));
      t(1) = getSlope(m2, m3, m4, m5);

      // right side: t(N-2), t(N-1)
      temp = (x(N - 2) - x(N - 1)) / (x(N - 3) - x(N - 1));
      temp = temp * temp;
      g0 = y(N - 1);
      g1 = ((y(N - 2) - y(N - 1)) - temp * (y(N - 3) - y(N - 1))) /
           ((x(N - 2) - x(N - 1)) - temp * (x(N - 3) - x(N - 1)));
      g2 = ((y(N - 3) - y(N - 1)) - g1 * (x(N - 3) - x(N - 1))) /
           ((x(N - 3) - x(N - 1)) * (x(N - 3) - x(N - 1)));
      x4 = x(N - 2) + (x(N - 1) - x(N - 3));
      x5 = x(N - 1) + (x(N - 1) - x(N - 3));
      y4 = g0 + g1 * (x4 - x(N - 1)) + g2 * (x4 - x(N - 1)) * (x4 - x(N - 1));
      y5 = g0 + g1 * (x5 - x(N - 1)) + g2 * (x5 - x(N - 1)) * (x5 - x(N - 1));
      m1 = (y(N - 3) - y(N - 4)) / (x(N - 3) - x(N - 4));
      m2 = (y(N - 2) - y(N - 3)) / (x(N - 2) - x(N - 3));
      m3 = (y(N - 1) - y(N - 2)) / (x(N - 1) - x(N - 2));
      m4 = (y4 - y(N - 1)) / (x4 - x(N - 1));
      m5 = (y5 - y4) / (x5 - x4);
      t(N - 2) = getSlope(m1, m2, m3, m4);
      t(N - 1) = getSlope(m2, m3, m4, m5);
      break;
    case splinePeriodic:
      // left: last two points determine the slopes t(0), t(1)
      m1 = (y(N - 1) - y(N - 2)) / (x(N - 1) - x(N - 2));
      m2 = (y(0) - y(N - 1)) / (x(0) - x(N - 1));
      m3 = (y(1) - y(0)) / (x(1) - x(0));
      m4 = (y(2) - y(1)) / (x(2) - x(1));
      m5 = (y(3) - y(2)) / (x(3) - x(2));
      t(0) = getSlope(m1, m2, m3, m4);
      t(1) = getSlope(m2, m3, m4, m5);
      // right: first two points determine the slopes t(N-2), t(N-1)
      m1 = (y(N - 3) - y(N - 4)) / (x(N - 3) - x(N - 4));
      m2 = (y(N - 2) - y(N - 3)) / (x(N - 2) - x(N - 3));
      m3 = (y(N - 1) - y(N - 2)) / (x(N - 1) - x(N - 2));
      m4 = (y(0) - y(N - 1)) / (x(0) - x(N - 1));
      m5 = (y(1) - y(0)) / (x(1) - x(0));
      t(N - 2) = getSlope(m1, m2, m3, m4);
      t(N - 1) = getSlope(m2, m3, m4, m5);
      break;
    case splineDerivativeZero:
      throw std::runtime_error(
          "erro in AkimaSpline::Interpolate: case splineDerivativeZero not "
          "implemented yet");
      break;
  }

  // calculate t's for all inner points [2,N-3]
  for (Index i = 2; i < N - 2; i++) {
    m1 = (y(i - 1) - y(i - 2)) / (x(i - 1) - x(i - 2));
    m2 = (y(i) - y(i - 1)) / (x(i) - x(i - 1));
    m3 = (y(i + 1) - y(i)) / (x(i + 1) - x(i));
    m4 = (y(i + 2) - y(i + 1)) / (x(i + 2) - x(i + 1));
    t(i) = getSlope(m1, m2, m3, m4);
  }

  // calculate p0,p1,p2,p3 for all intervals 0..(N-2), where interval
  // [x(i),x(i+1)] shall have number i (this means that the last interval
  // has number N-2)
  for (Index i = 0; i < N - 1; i++) {
    p0(i) = y(i);
    p1(i) = t(i);
    p2(i) =
        (3.0 * (y(i + 1) - y(i)) / (x(i + 1) - x(i)) - 2.0 * t(i) - t(i + 1)) /
        (x(i + 1) - x(i));
    p3(i) = (t(i) + t(i + 1) - 2.0 * (y(i + 1) - y(i)) / (x(i + 1) - x(i))) /
            ((x(i + 1) - x(i)) * (x(i + 1) - x(i)));
  }
}

void AkimaSpline::Fit(Eigen::VectorXd &, Eigen::VectorXd &) {
  throw std::runtime_error("Akima fit not implemented.");
}

}  // namespace tools
}  // namespace votca
