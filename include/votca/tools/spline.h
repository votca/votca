/* 
 * File:   spline.h
 * Author: hahng
 *
 * Created on December 10, 2010, 1:42 PM
 */

#ifndef __VOTCA_SPLINE_H
#define	__VOTCA_SPLINE_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>

namespace votca{namespace tools{

namespace ub = boost::numeric::ublas;

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
    virtual void Interpolate(ub::vector<double> &x, ub::vector<double> &y) = 0;
    
    /**
     * \brief Fit spline through noisy (x,y) values. Points on resulting fitted spline can be obtained via Calculate().
     * \param x values of data to be fitted
     * \param y values of data to be fitted
     * both vectors must be of same size
     */
    virtual void Fit(ub::vector<double> &x, ub::vector<double> &y) = 0;

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
        splinePeriodic,    ///< periodic boundary conditions: \f$f_0=f_N\f$
        splineDerivativeZero ///< derivatives and end-points are zero.
    };

    /**
     * \brief Set the boundary type of the spline
     * \param boundary of type eBoundary
     */
    void setBC(eBoundary bc) {_boundaries = bc;}

    /**
     * \brief Get the grid point of certain index
     * \param index of grid point
     * \return grid value
     */
    inline double getGridPoint(const size_t &i);

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
    ub::vector<double> &getX() {return _r; }

    /**
     * \brief Get the spline data _f
     * \return pointer to the corresponding array
     */
    ub::vector<double> &getSplineF() { return _f; }

    /**
     * \brief Get second derivatives (cubic splines)
     * \return pointer to the corresponding array
     */
    ub::vector<double> &getSplineF2() { return _f2; }
    
protected:
    eBoundary _boundaries;
    // the grid points
    ub::vector<double> _r;
    // y values of grid points
    ub::vector<double> _f;
    // second derivatives of grid points
    ub::vector<double> _f2;
};

template<typename vector_type1, typename vector_type2>
inline void Spline::Calculate(vector_type1 &x, vector_type2 &y)
{
    for(size_t i=0; i<x.size(); ++i)
        y(i) = Calculate(x(i));
}

template<typename vector_type1, typename vector_type2>
inline void Spline::CalculateDerivative(vector_type1 &x, vector_type2 &y)
{
    for(size_t i=0; i<x.size(); ++i)
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
    size_t i;
    for(i=0; i<_r.size(); ++i)
        if(_r[i]>r) break;
    return i-1;
}

inline double Spline::getGridPoint(const size_t &i)
{
    if(i>=_r.size()) {
        return 0;
    }
    return _r[i];
}

}}

#endif	/* __VOTCA_SPLINE_H */
