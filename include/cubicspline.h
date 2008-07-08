/* 
 * File:   CubicSpline.h
 * Author: ruehle
 *
 * Created on July 4, 2008, 3:03 PM
 */

#ifndef _CUBICSPLINE_H
#define	_CUBICSPLINE_H

#include <boost/numeric/ublas/vector.hpp>
namespace ub = boost::numeric::ublas;

class CubicSpline
{
public:
    CubicSpline() {}
    ~CubicSpline() {}
    
    // Generates the r_k
    void GenerateGrid(double &min, double &max, double &h);
    
    double A(double r);
    double B(double r);
    double C(double r);
    double D(double r);
    
protected:
    
    // in which interval is point r, return i for interval r_i r_{i+1}, -1 for out of range
    int getInterval(double &r);
    ub::vector<double> _r;
};

inline void CubicSpline::GenerateGrid(double &min, double &max, double &h)
{
    for (double r_init = min; r_init < max; r_init += h ) {
            _r.push_back(r_init);
    }
}

inline double CubicSpline::A(double &r)
{
    return ( 1.0 - (r -_r[getInverval(r)])/(_r[getInterval(r)+1]-_r[getInterval(r)]) );
}

inline double CubicSpline::B(double r)
{
    return  (r -_r[getInverval(r)])/(_r[getInterval(r)+1]-_r[getInterval(r)]) ;
}

inline double CubicSpline::C(double r)
{
    double xxi, h;
    xxi = r -_r[getInverval(r)];
    h   = _r[getInterval(r)+1]-_r[getInterval(r)];
    
    return ( 0.5*xxi*xxi - (1.0/6.0)*xxi*xxi*xxi/h - (1.0/3.0)*xxi*h) ;
}

inline double CubicSpline::D(double r)
{
    double xxi, h;
    xxi = r -_r[getInverval(r)];
    h   = _r[getInterval(r)+1]-_r[getInterval(r)]; 
    
    return ( (1.0/6.0)*xxi*xxi*xxi/h - (1.0/6.0)*xxi*h ) ;
}

inline int CubicSpline::getInterval(double &r)
{
    if (r > _r.back || r < _r.front) return -1;
    return int( (r - _r.front()) / (_r.back() - _r.front()) * _r.size() );
}


#endif	/* _CUBICSPLINE_H */

