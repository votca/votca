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

#ifndef VOTCA_TOOLS_CUBICSPLINE_H
#define VOTCA_TOOLS_CUBICSPLINE_H

#include "spline.h"
#include <iostream>
#include <votca/tools/eigen.h>

namespace votca {
namespace tools {

/**
    \brief A cubic spline class

    This class does cubic piecewise spline interpolation and spline fitting.
    As representation of a single spline, the general form
    \f[
        S_i(x) = A(x,h_i) f_i + B(x,h_i) f_{i+1} + C(x,h_i) f''_i + d(x,h_i)
   f''_{i+1} \f] with \f[
        x_i \le x < x_{i+1}\,,\\
        h_i = x_{i+1} - x_{i}
    \f]
    The \f$f_i\,,\,,f''_i\f$ are the function values and second derivates
    at point \f$x_i\f$.

    The parameters \f$f''_i\f$ are no free parameters, they are determined by
   the smoothing condition that the first derivatives are continuous. So the
   only free paramers are the grid points x_i as well as the functon values f_i
   at these points. A spline can be generated in several ways:
    - Interpolation spline
    - Fitting spline (fit to noisy data)
    - calculate the parameters elsewere and fill the spline class
*/

class CubicSpline : public Spline {
 public:
  // default constructor
  CubicSpline() = default;
  // CubicSpline() :
  //    _boundaries(splineNormal) {}

  // destructor
  ~CubicSpline() override = default;

  // construct an interpolation spline
  // x, y are the the points to construct interpolation, both vectors must be of
  // same size
  void Interpolate(Eigen::VectorXd &x, Eigen::VectorXd &y) override;

  // fit spline through noisy data
  // x,y are arrays with noisy data, both vectors must be of same size
  void Fit(Eigen::VectorXd &x, Eigen::VectorXd &y) override;

  // Calculate the function value
  double Calculate(const double &x) override;

  // Calculate the function derivative
  double CalculateDerivative(const double &x) override;

  // Calculate the function value for a whole array, story it in y
  template <typename vector_type1, typename vector_type2>
  void Calculate(vector_type1 &x, vector_type2 &y);

  // Calculate the derivative value for a whole array, story it in y
  template <typename vector_type1, typename vector_type2>
  void CalculateDerivative(vector_type1 &x, vector_type2 &y);

  // set spline parameters to values that were externally computed
  template <typename vector_type>
  void setSplineData(vector_type &f, vector_type &f2) {
    _f = f;
    _f2 = f2;
  }

  /**
   * \brief Add a point (one entry) to fitting matrix
   * \param pointer to matrix
   * \param value x
   * \param offsets relative to getInterval(x)
   * \param scale parameters for terms "A,B,C,D"
   * When creating a matrix to fit data with a spline, this function creates
   * one entry in that fitting matrix.
   */
  template <typename matrix_type>
  void AddToFitMatrix(matrix_type &A, double x, long offset1, long offset2 = 0,
                      double scale = 1);

  /**
   * \brief Add a point (one entry) to fitting matrix
   * \param pointer to matrix [in] [out]
   * \param value x [in]
   * \param offsets relative to getInterval(x) [in]
   * \param scale1 parameters for terms "A,B,C,D" [in]
   * \param scale2 parameters for terms "AA,BB,CC,DD" [in]
   * When creating a matrix to fit data with a spline, this function creates
   * one entry in that fitting matrix.
   */
  template <typename matrix_type>
  void AddToFitMatrix(matrix_type &A, double x, int long offset1,
                      int long offset2, double scale1, double scale2);

  /**
   * \brief Add a vector of points to fitting matrix
   * \param pointer to matrix
   * \param vector of x values
   * \param offsets relative to getInterval(x)
   * Same as previous function, but vector-valued and with scale=1.0
   */
  template <typename matrix_type, typename vector_type>
  void AddToFitMatrix(matrix_type &M, vector_type &x, long offset1,
                      long offset2 = 0);

  /**
   * \brief Add boundary condition of sum_i f_i =0 to fitting matrix
   * \param pointer to matrix
   * \param offsets
   */
  template <typename matrix_type>
  void AddBCSumZeroToFitMatrix(matrix_type &A, long offset1, long offset2 = 0);

  /**
   * \brief Add boundary conditions to fitting matrix
   * \param pointer to matrix
   * \param offsets
   */
  template <typename matrix_type>
  void AddBCToFitMatrix(matrix_type &A, long offset1, long offset2 = 0);

 protected:
  // A spline can be written in the form
  // S_i(x) =   A(x,x_i,x_i+1)*f_i     + B(x,x_i,x_i+1)*f''_i
  //          + C(x,x_i,x_i+1)*f_{i+1} + D(x,x_i,x_i+1)*f''_{i+1}
  double A(const double &r);
  double B(const double &r);
  double C(const double &r);
  double D(const double &r);

  double Aprime(const double &r);
  double Bprime(const double &r);
  double Cprime(const double &r);
  double Dprime(const double &r);

  // tabulated derivatives at grid points. Second argument: 0 - left, 1 - right
  double A_prime_l(long i);
  double A_prime_r(long i);
  double B_prime_l(long i);
  double B_prime_r(long i);
  double C_prime_l(long i);
  double C_prime_r(long i);
  double D_prime_l(long i);
  double D_prime_r(long i);
};

inline double CubicSpline::Calculate(const double &r) {
  long interval = getInterval(r);
  return A(r) * _f[interval] + B(r) * _f[interval + 1] + C(r) * _f2[interval] +
         D(r) * _f2[interval + 1];
}

inline double CubicSpline::CalculateDerivative(const double &r) {
  long interval = getInterval(r);
  return Aprime(r) * _f[interval] + Bprime(r) * _f[interval + 1] +
         Cprime(r) * _f2[interval] + Dprime(r) * _f2[interval + 1];
}

template <typename matrix_type>
inline void CubicSpline::AddToFitMatrix(matrix_type &M, double x, long offset1,
                                        long offset2, double scale) {
  long spi = getInterval(x);
  M(offset1, offset2 + spi) += A(x) * scale;
  M(offset1, offset2 + spi + 1) += B(x) * scale;
  M(offset1, offset2 + spi + _r.size()) += C(x) * scale;
  M(offset1, offset2 + spi + _r.size() + 1) += D(x) * scale;
}

// for adding f'(x)*scale1 + f(x)*scale2 as needed for threebody interactions
template <typename matrix_type>
inline void CubicSpline::AddToFitMatrix(matrix_type &M, double x, long offset1,
                                        long offset2, double scale1,
                                        double scale2) {
  long spi = getInterval(x);
  M(offset1, offset2 + spi) += Aprime(x) * scale1;
  M(offset1, offset2 + spi + 1) += Bprime(x) * scale1;
  M(offset1, offset2 + spi + _r.size()) += Cprime(x) * scale1;
  M(offset1, offset2 + spi + _r.size() + 1) += Dprime(x) * scale1;

  M(offset1, offset2 + spi) += A(x) * scale2;
  M(offset1, offset2 + spi + 1) += B(x) * scale2;
  M(offset1, offset2 + spi + _r.size()) += C(x) * scale2;
  M(offset1, offset2 + spi + _r.size() + 1) += D(x) * scale2;
}

template <typename matrix_type, typename vector_type>
inline void CubicSpline::AddToFitMatrix(matrix_type &M, vector_type &x,
                                        long offset1, long offset2) {
  for (int i = 0; i < x.size(); ++i) {
    long spi = getInterval(x(i));
    M(offset1 + i, offset2 + spi) = A(x(i));
    M(offset1 + i, offset2 + spi + 1) = B(x(i));
    M(offset1 + i, offset2 + spi + _r.size()) = C(x(i));
    M(offset1 + i, offset2 + spi + _r.size() + 1) = D(x(i));
  }
}

template <typename matrix_type>
inline void CubicSpline::AddBCSumZeroToFitMatrix(matrix_type &M, long offset1,
                                                 long offset2) {
  for (int i = 0; i < _r.size(); ++i) {
    M(offset1, offset2 + i) = 1;
    M(offset1, offset2 + _r.size() + i) = 0;
  }
}

template <typename matrix_type>
inline void CubicSpline::AddBCToFitMatrix(matrix_type &M, long offset1,
                                          long offset2) {
  for (int i = 0; i < _r.size() - 2; ++i) {
    M(offset1 + i + 1, offset2 + i) = A_prime_l(i);
    M(offset1 + i + 1, offset2 + i + 1) = B_prime_l(i) - A_prime_r(i);
    M(offset1 + i + 1, offset2 + i + 2) = -B_prime_r(i);

    M(offset1 + i + 1, offset2 + _r.size() + i) = C_prime_l(i);
    M(offset1 + i + 1, offset2 + _r.size() + i + 1) =
        D_prime_l(i) - C_prime_r(i);
    M(offset1 + i + 1, offset2 + _r.size() + i + 2) = -D_prime_r(i);
  }
  // currently only natural boundary conditions:
  switch (_boundaries) {
    case splineNormal:
      M(offset1, offset2 + _r.size()) = 1;
      M(offset1 + _r.size() - 1, offset2 + 2 * _r.size() - 1) = 1;
      break;
    case splineDerivativeZero:
      // y
      M(offset1 + 0, offset2 + 0) = -1 * A_prime_l(0);
      M(offset1 + 0, offset2 + 1) = -1 * B_prime_l(0);

      M(offset1 + _r.size() - 1, offset2 + _r.size() - 2) =
          A_prime_l(_r.size() - 2);
      M(offset1 + _r.size() - 1, offset2 + _r.size() - 1) =
          B_prime_l(_r.size() - 2);

      // y''
      M(offset1 + 0, offset2 + _r.size() + 0) = D_prime_l(0);
      M(offset1 + 0, offset2 + _r.size() + 1) = C_prime_l(0);

      M(offset1 + _r.size() - 1, offset2 + 2 * _r.size() - 2) =
          C_prime_l(_r.size() - 2);
      M(offset1 + _r.size() - 1, offset2 + 2 * _r.size() - 1) =
          D_prime_l(_r.size() - 2);
      break;

    case splinePeriodic:
      M(offset1, offset2) = 1;
      M(offset1, offset2 + _r.size() - 1) = -1;
      M(offset1 + _r.size() - 1, offset2 + _r.size()) = 1;
      M(offset1 + _r.size() - 1, offset2 + 2 * _r.size() - 1) = -1;
      break;
  }
}

inline double CubicSpline::A(const double &r) {
  return (1.0 - (r - _r[getInterval(r)]) /
                    (_r[getInterval(r) + 1] - _r[getInterval(r)]));
}

inline double CubicSpline::Aprime(const double &r) {
  return -1.0 / (_r[getInterval(r) + 1] - _r[getInterval(r)]);
}

inline double CubicSpline::B(const double &r) {
  return (r - _r[getInterval(r)]) /
         (_r[getInterval(r) + 1] - _r[getInterval(r)]);
}

inline double CubicSpline::Bprime(const double &r) {
  return 1.0 / (_r[getInterval(r) + 1] - _r[getInterval(r)]);
}

inline double CubicSpline::C(const double &r) {
  double xxi, h;
  xxi = r - _r[getInterval(r)];
  h = _r[getInterval(r) + 1] - _r[getInterval(r)];

  return (0.5 * xxi * xxi - (1.0 / 6.0) * xxi * xxi * xxi / h -
          (1.0 / 3.0) * xxi * h);
}

inline double CubicSpline::Cprime(const double &r) {
  double xxi, h;
  xxi = r - _r[getInterval(r)];
  h = _r[getInterval(r) + 1] - _r[getInterval(r)];

  return (xxi - 0.5 * xxi * xxi / h - h / 3);
}

inline double CubicSpline::D(const double &r) {
  double xxi, h;
  xxi = r - _r[getInterval(r)];
  h = _r[getInterval(r) + 1] - _r[getInterval(r)];

  return ((1.0 / 6.0) * xxi * xxi * xxi / h - (1.0 / 6.0) * xxi * h);
}

inline double CubicSpline::Dprime(const double &r) {
  double xxi, h;
  xxi = r - _r[getInterval(r)];
  h = _r[getInterval(r) + 1] - _r[getInterval(r)];

  return (0.5 * xxi * xxi / h - (1.0 / 6.0) * h);
}

/**
inline int CubicSpline::getInterval(double &r)
{
    if (r < _r[0] || r > _r[_r.size() - 1]) return -1;
    return int( (r - _r[0]) / (_r[_r.size()-1] - _r[0]) * (_r.size() - 1) );
}
 **/

inline double CubicSpline::A_prime_l(long i) {
  return -1.0 / (_r[i + 1] - _r[i]);
}

inline double CubicSpline::B_prime_l(long i) {
  return 1.0 / (_r[i + 1] - _r[i]);
}

inline double CubicSpline::C_prime_l(long i) {
  return (1.0 / 6.0) * (_r[i + 1] - _r[i]);
}

inline double CubicSpline::D_prime_l(long i) {
  return (1.0 / 3.0) * (_r[i + 1] - _r[i]);
}

inline double CubicSpline::A_prime_r(long i) {
  return -1.0 / (_r[i + 2] - _r[i + 1]);
}

inline double CubicSpline::B_prime_r(long i) {
  return 1.0 / (_r[i + 2] - _r[i + 1]);
}

inline double CubicSpline::C_prime_r(long i) {
  return -(1.0 / 3.0) * (_r[i + 2] - _r[i + 1]);
}

inline double CubicSpline::D_prime_r(long i) {
  return -(1.0 / 6.0) * (_r[i + 2] - _r[i + 1]);
}

}  // namespace tools
}  // namespace votca

#endif /* VOTCA_TOOLS_CUBICSPLINE_H */
