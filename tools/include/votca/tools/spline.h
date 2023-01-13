/*
 * Copyright 2010-2020 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_TOOLS_SPLINE_H
#define VOTCA_TOOLS_SPLINE_H

// Local VOTCA includes
#include "eigen.h"
#include "types.h"

namespace votca {
namespace tools {

/**
 * \brief Spline Class
 *
 *  class supports spline interpolation and fit of data
 *  the cubic spline class, akima spline class and linear spline class are
 * inherited from this one
 */

class Spline {
 public:
  Spline() = default;

  virtual ~Spline() = default;

  /**
   * \brief Calculate interpolating spline for given (x,y) values. Points on
   * resulting spline can be obtained via Calculate(). \param x values of data
   * to be interpolated \param y values of data to be interpolated both vectors
   * must be of same size
   */
  virtual void Interpolate(const Eigen::VectorXd &x,
                           const Eigen::VectorXd &y) = 0;

  /**
   * \brief Fit spline through noisy (x,y) values. Points on resulting fitted
   * spline can be obtained via Calculate(). \param x values of data to be
   * fitted \param y values of data to be fitted both vectors must be of same
   * size
   */
  virtual void Fit(const Eigen::VectorXd &x, const Eigen::VectorXd &y) = 0;

  /**
   * \brief Calculate spline function value for a given x value on the spline
   * created by Interpolate() or Fit() \param x data value \return y value
   */
  virtual double Calculate(double x) = 0;

  /**
   * \brief Calculate y value for a given x value on the derivative of the
   * spline created by function Interpolate or Fit \param x data value \return y
   * value of derivative
   */
  virtual double CalculateDerivative(double x) = 0;

  /// enum for type of boundary condition
  enum eBoundary {
    splineNormal = 0,         ///< normal boundary conditions: \f$f_0=f_N=0\f$
    splinePeriodic = 1,       ///< periodic boundary conditions: \f$f_0=f_N\f$
    splineDerivativeZero = 2  ///< derivatives and end-points are zero.
  };

  /**
   * \brief Set the boundary type of the spline
   * \param bc of type eBoundary
   */
  void setBC(eBoundary bc) { boundaries_ = bc; }

  /**
   * \brief Set the boundary type of the spline
   * \param bc of type Index
   */
  void setBCInt(Index bc) {
    switch (bc) {
      case 0:
        boundaries_ = splineNormal;
        break;
      case 1:
        boundaries_ = splinePeriodic;
        break;
      case 2:
        boundaries_ = splineDerivativeZero;
        break;
    }
  }

  /**
   * \brief Get the grid point of certain index
   * \param i index of grid point
   * \return grid value
   */
  double getGridPoint(int i);

  /**
   * \brief Calculate spline function values for given x values on the spline
   * created by Interpolate() or Fit()
   * \param x vector of data values
   * \return vector of y value
   */
  Eigen::VectorXd Calculate(const Eigen::VectorXd &x);

  /**
   * \brief Calculate y values for given x values on the derivative of the
   * spline created by function Interpolate or Fit
   * \param x vector of data values
   * \return vector of y value
   */
  Eigen::VectorXd CalculateDerivative(const Eigen::VectorXd &x);

  /**
   * \brief Print spline values (using Calculate()) on output "out" on the
   * entire grid in steps of "interval"
   * \param out output
   * \param interval size of "interval"
   */
  void Print(std::ostream &out, double interval);

  /**
   * \brief Determine the index of the interval containing value r
   * \param r
   * \return interval index
   */
  Index getInterval(double r);

  /**
   * \brief Generate the grid for fitting from "min" to "max" in steps of "h"
   * \param min left interval border
   * \param max right interval border
   * \param h step
   * \return number of grid values in the interval
   */
  Index GenerateGrid(double min, double max, double h);

  /**
   * \brief Get the grid array x
   * \return pointer to the corresponding array
   */
  Eigen::VectorXd &getX() { return r_; }
  const Eigen::VectorXd &getX() const { return r_; }
  /**
   * \brief Get the spline data  f_
   * \return reference to the corresponding array
   */
  // Eigen::VectorXd &getSplineF() { return  f_; }
  // const Eigen::VectorXd &getSplineF() const { return  f_; }

  /**
   * \brief Get second derivatives (cubic splines)
   * \return reference to the corresponding array
   */
  // Eigen::VectorXd &getSplineF2() { return  f2_; }
  // const Eigen::VectorXd &getSplineF2() const { return  f_; }

 protected:
  eBoundary boundaries_ = eBoundary::splineNormal;
  // the grid points
  Eigen::VectorXd r_;
};

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_SPLINE_H
