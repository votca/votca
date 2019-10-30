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

#include <cmath>
#include <iostream>
#include <votca/tools/cubicspline.h>
#include <votca/tools/linalg.h>

namespace votca {
namespace tools {

using namespace std;

void CubicSpline::Interpolate(const Eigen::VectorXd &x,
                              const Eigen::VectorXd &y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument(
        "error in CubicSpline::Interpolate : sizes of vectors x and y do not "
        "match");
  }

  if (x.size() < 3) {
    throw std::invalid_argument(
        "error in CubicSpline::Interpolate : vectors x and y have to contain "
        "at least 3 points");
  }

  const Index N = x.size();

  // copy the grid points into f
  _r = x;
  _f = y;
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(N);

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);

  for (Index i = 0; i < N - 2; ++i) {
    temp(i + 1) =
        -(A_prime_l(i) * _f(i) + (B_prime_l(i) - A_prime_r(i)) * _f(i + 1) -
          B_prime_r(i) * _f(i + 2));

    A(i + 1, i) = C_prime_l(i);
    A(i + 1, i + 1) = D_prime_l(i) - C_prime_r(i);
    A(i + 1, i + 2) = -D_prime_r(i);
  }

  switch (_boundaries) {
    case splineNormal:
      A(0, 0) = 1;
      A(N - 1, N - 1) = 1;
      break;
    case splinePeriodic:
      A(0, 0) = 1;
      A(0, N - 1) = -1;
      A(N - 1, 0) = 1;
      A(N - 1, N - 1) = -1;
      break;
    case splineDerivativeZero:
      throw std::runtime_error(
          "erro in CubicSpline::Interpolate: case splineDerivativeZero not "
          "implemented yet");
      break;
  }

  Eigen::HouseholderQR<Eigen::MatrixXd> QR(A);
  _f2 = QR.solve(temp);
}

void CubicSpline::Fit(const Eigen::VectorXd &x, const Eigen::VectorXd &y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument(
        "error in CubicSpline::Fit : sizes of vectors x and y do not match");
  }

  const Index N = x.size();
  const Index ngrid = _r.size();

  // construct the equation
  // A*u = b
  // where u = { {f[i]}, {f''[i]} }
  // and b[i] = y[i] for 0<=i<N
  // and b[i]=0 for i>=N (for smoothing condition)
  // A[i,j] contains the data fitting + the spline smoothing conditions

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, 2 * ngrid);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(
      ngrid, 2 * ngrid);  // Matrix with smoothing conditions

  // Construct smoothing matrix
  AddBCToFitMatrix(B, 0);
  // construct the matrix to fit the points and the vector b
  AddToFitMatrix(A, x, 0);
  // now do a constrained qr solve
  Eigen::VectorXd sol = linalg_constrained_qrsolve(A, y, B);

  // check vector "sol" for nan's
  for (Index i = 0; i < 2 * ngrid; i++) {
    if ((std::isinf(sol(i))) || (std::isnan(sol(i)))) {
      throw std::runtime_error(
          "error in CubicSpline::Fit : value nan occurred due to wrong fitgrid "
          "boundaries");
    }
  }

  _f = sol.segment(0, ngrid);
  _f2 = sol.segment(ngrid, ngrid);
}

}  // namespace tools
}  // namespace votca
