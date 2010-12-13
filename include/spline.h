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

class Spline
{
public:
    Spline() :
        _boundaries(splineNormal) {}

    virtual ~Spline() {}

    /**
     * \brief executes the program
     * \param argc argc from main
     * \param argv argv from main
     * \return return code
     */
    virtual void Interpolate(ub::vector<double> &x, ub::vector<double> &y) = 0;

    virtual void Fit(ub::vector<double> &x, ub::vector<double> &y) = 0;
    virtual double Calculate(const double &x) = 0;
    virtual double CalculateDerivative(const double &x) = 0;

    /// enum for type of boundary condition
    enum eBoundary {
        splineNormal = 0,  ///< normal boundary conditions: \f$f_0=f_N=0\f$
        splinePeriodic,    ///< periodic boundary conditions: \f$f_0=f_N\f$
        splineDerivativeZero ///< derivatives and end-points are zero.
    };
    /// set the boundary type of the spline
    void setBC(eBoundary bc) {_boundaries = bc;}

    /// ERROR-PRONE implementation, make it better!!!
    double getGridPoint(const int &i) {return _r[i];}

    template<typename vector_type1, typename vector_type2>
    inline void Calculate(vector_type1 &x, vector_type2 &y);

    template<typename vector_type1, typename vector_type2>
    inline void CalculateDerivative(vector_type1 &x, vector_type2 &y);

    inline void Print(std::ostream &out, double interval);

    inline int getInterval(const double &r);

    int GenerateGrid(double min, double max, double h);

    
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

}}

#endif	/* __VOTCA_SPLINE_H */