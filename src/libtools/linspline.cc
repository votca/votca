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
#include <votca/tools/linalg.h>
#include <votca/tools/linspline.h>

namespace votca {
namespace tools {

using namespace std;

void LinSpline::Interpolate(Eigen::VectorXd &x, Eigen::VectorXd &y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument(
        "error in LinSpline::Interpolate : sizes of vectors x and y do not "
        "match");
  }

  if (x.size() < 2) {
    throw std::invalid_argument(
        "error in LinSpline::Interpolate : vectors x and y have to contain at "
        "least 2 points");
  }

  const long int N = x.size();

  // adjust the grid
  _r.resize(N);

  // copy the grid points into f
  _r = x;

  // LINEAR SPLINE: a(i) * x + b(i)
  // where i=number of interval

  // initialize vectors a,b
  a = Eigen::VectorXd::Zero(N);
  b = Eigen::VectorXd::Zero(N);

  // boundary conditions not applicable

  // calculate a,b for all intervals 0..(N-2), where interval
  // [x(i),x(i+1)] shall have number i (this means that the last interval
  // has number N-2)
  for (int i = 0; i < N - 1; i++) {
    a(i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));
    b(i) = y(i) - a(i) * x(i);
  }
}

void LinSpline::Fit(Eigen::VectorXd &x, Eigen::VectorXd &y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument(
        "error in LinSpline::Fit : sizes of vectors x and y do not match");
  }

  const long int N = x.size();
  const long int ngrid = _r.size();

  // construct the equation
  // A*u = b
  // The matrix A contains all conditions
  // s_i(x) = (y(i+1)-y(i)) * (x-r(i))/(r(i+1)-r(i)) + y(i)
  // where y(i) are the unknown values at grid points r(i), and
  // the condition y=s_i(x) is to be satisfied at all input points:
  // therefore b=y and u=vector of all unknown y(i)

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, ngrid);

  // construct matrix A
  for (int i = 0; i < N; i++) {
    long int interval = getInterval(x(i));
    A(i, interval) =
        1 - (x(i) - _r(interval)) / (_r(interval + 1) - _r(interval));
    A(i, interval + 1) =
        (x(i) - _r(interval)) / (_r(interval + 1) - _r(interval));
  }

  // now do a qr solve
  Eigen::HouseholderQR<Eigen::MatrixXd> QR(A);
  Eigen::VectorXd sol = QR.solve(y);

  // vector "sol" contains all y-values of fitted linear splines at each
  // interval border
  // get a(i) and b(i) for piecewise splines out of solution vector "sol"
  a = Eigen::VectorXd::Zero(ngrid - 1);
  b = Eigen::VectorXd::Zero(ngrid - 1);
  for (int i = 0; i < ngrid - 1; i++) {
    a(i) = (sol(i + 1) - sol(i)) / (_r(i + 1) - _r(i));
    b(i) = -a(i) * _r(i) + sol(i);
  }
}
}  // namespace tools
}  // namespace votca
