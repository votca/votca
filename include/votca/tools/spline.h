/* 
 * Copyright 2010-2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_SPLINE_H
#define	__VOTCA_SPLINE_H

#include <votca/tools/eigen.h>

namespace votca{namespace tools{


/**
 * \brief Spline Class
 *
 *  class supports spline interpolation and fit of data
 *  the cubic spline class, akima spline class and linear spline class are inherited from this one
*/

class Spline
{
public:
    Spline() :
        _boundaries(splineNormal) {}

    virtual ~Spline() {}

    /**
     * \brief Calculate interpolating spline for given (x,y) values. Points on resulting spline can be obtained via Calculate().
     * \param x values of data to be interpolated
     * \param y values of data to be interpolated
     * both vectors must be of same size
     */
    virtual void Interpolate(Eigen::VectorXd &x, Eigen::VectorXd &y) = 0;
    
    /**
     * \brief Fit spline through noisy (x,y) values. Points on resulting fitted spline can be obtained via Calculate().
     * \param x values of data to be fitted
     * \param y values of data to be fitted
     * both vectors must be of same size
     */
    virtual void Fit(Eigen::VectorXd &x, Eigen::VectorXd &y) = 0;

    /**
     * \brief Calculate spline function value for a given x value on the spline created by Interpolate() or Fit()
     * \param x data value
     * \return y value
     */
    virtual double Calculate(const double &x) = 0;

    /**
     * \brief Calculate y value for a given x value on the derivative of the spline created by function Interpolate or Fit
     * \param x data value
     * \return y value of derivative
     */
    virtual double CalculateDerivative(const double &x) = 0;

    /// enum for type of boundary condition
    enum eBoundary {
        splineNormal = 0,  ///< normal boundary conditions: \f$f_0=f_N=0\f$
        splinePeriodic = 1,    ///< periodic boundary conditions: \f$f_0=f_N\f$
        splineDerivativeZero = 2 ///< derivatives and end-points are zero.
    };

    /**
     * \brief Set the boundary type of the spline
     * \param boundary of type eBoundary
     */
    void setBC(eBoundary bc) {_boundaries = bc;}
    
        /**
     * \brief Set the boundary type of the spline
     * \param boundary of type int
     */
    void setBCInt(int bc) {        
        switch(bc){
            case 0:
                _boundaries = splineNormal;
                break;
            case 1:
                _boundaries = splinePeriodic;
                break;
            case 2:
                _boundaries = splineDerivativeZero;
                break;                    
        }         
    }

    /**
     * \brief Get the grid point of certain index
     * \param index of grid point
     * \return grid value
     */
    inline double getGridPoint(const int &i);

    /**
     * \brief Calculate spline function values for given x values on the spline created by Interpolate() or Fit()
     * \param vector of x data values
     * \return vector of y value
     */
    template<typename vector_type1, typename vector_type2>
    inline void Calculate(vector_type1 &x, vector_type2 &y);

    /**
     * \brief Calculate y values for given x values on the derivative of the spline created by function Interpolate or Fit
     * \param vector of x data values
     * \return vector of y value
     */
    template<typename vector_type1, typename vector_type2>
    inline void CalculateDerivative(vector_type1 &x, vector_type2 &y);

    /**
     * \brief Print spline values (using Calculate()) on output "out" on the entire grid in steps of "interval"
     * \param reference "out" to output
     * \param steps of size "interval"
     */
    inline void Print(std::ostream &out, double interval);

    /**
     * \brief Determine the index of the interval containing value r
     * \param value r
     * \return interval index
     */
    inline int getInterval(const double &r);

    /**
     * \brief Generate the grid for fitting from "min" to "max" in steps of "h"
     * \param left interval border "min"
     * \param right interval border "max"
     * \param step "h"
     * \return number of grid values in the interval
     */
    int GenerateGrid(double min, double max, double h);

    /**
     * \brief Get the grid array x
     * \return pointer to the corresponding array
     */
    Eigen::VectorXd &getX() {return _r; }

    /**
     * \brief Get the spline data _f
     * \return pointer to the corresponding array
     */
    Eigen::VectorXd &getSplineF() { return _f; }

    /**
     * \brief Get second derivatives (cubic splines)
     * \return pointer to the corresponding array
     */
    Eigen::VectorXd &getSplineF2() { return _f2; }
    
protected:
    eBoundary _boundaries;
    // the grid points
    Eigen::VectorXd _r;
    // y values of grid points
    Eigen::VectorXd _f;
    // second derivatives of grid points
    Eigen::VectorXd _f2;
};

template<typename vector_type1, typename vector_type2>
inline void Spline::Calculate(vector_type1 &x, vector_type2 &y)
{
    for(int i=0; i<x.size(); ++i)
        y(i) = Calculate(x(i));
}

template<typename vector_type1, typename vector_type2>
inline void Spline::CalculateDerivative(vector_type1 &x, vector_type2 &y)
{
    for(int i=0; i<x.size(); ++i)
        y(i) = CalculateDerivative(x(i));
}

inline void Spline::Print(std::ostream &out, double interval)
{
    for (double x = _r[0]; x < _r[_r.size() - 1]; x += interval)
        out << x << " " << Calculate(x) << "\n";
}

inline int Spline::getInterval(const double &r)
{
    if (r < _r[0]) return 0;
    if(r > _r[_r.size() - 2]) return _r.size()-2;
    int i;
    for( i=0; i<_r.size(); ++i)
        if(_r[i]>r) break;
    return i-1;
}

inline double Spline::getGridPoint(const int &i)
{
    if(i>=_r.size()) {
        return 0;
    }
    return _r[i];
}

}}

#endif	/* __VOTCA_SPLINE_H */
