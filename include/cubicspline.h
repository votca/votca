/* 
 * File:   CubicSpline.h
 * Author: ruehle
 *
 * Created on July 4, 2008, 3:03 PM
 */

#ifndef _CUBICSPLINE_H
#define	_CUBICSPLINE_H

#include <boost/numeric/ublas/vector.hpp>
#include <iostream>

namespace ub = boost::numeric::ublas;

/**
    \brief A cubic spline class
  
    This class represents a Cubic spline. Spline interpolation
    as well as spline fitting can be done.
 */

class CubicSpline
{
public:
    // default constructor
    CubicSpline() :
        _boundaries(splineNormal) {}
    
    // destructor
    ~CubicSpline() {};    
    
    // enum for type of boundary condition
    enum eBoundary {
        splineNormal = 0,
        splinePeriodic
    };
    
    //void Interpolation(ub::vector<double> x, ub::vector<double> y);
    
    void Fit(ub::vector<double> x, ub::vector<double> y);
    
    // Generates the r_k, returns the number of grid points
    int GenerateGrid(double min, double max, double h);
    

    // give string in rangeparser format: e.g. 1:0.1:10;11:1:20
    //int GenerateGrid(string range);
    
    // A spline can be written in the form
    // S_i(x) =   A(x,x_i,x_i+1)*f_i     + B(x,x_i,x_i+1)*f''_i 
    //          + C(x,x_i,x_i+1)*f_{i+1} + D(x,x_i,x_i+1)*f''_{i+1}
    double A(double &r);
    double B(double &r);
    double C(double &r);
    double D(double &r);
  
    // tabulated derivatives at grid points. Second argument: 0 - left, 1 - right
    double A_prime(int i, int flag); 
    double B_prime(int i, int flag);    
    double C_prime(int i, int flag);
    double D_prime(int i, int flag);
    
    double A_prime_l(int i);     
    double A_prime_r(int i);     
    double B_prime_l(int i);    
    double B_prime_r(int i);    
    double C_prime_l(int i);
    double C_prime_r(int i);
    double D_prime_l(int i);
    double D_prime_r(int i);
    

    // store data in the spline, todo: rename that function
    void GetResult(ub::vector<double> *x_pointer);
    double getFunctionValue(double &r);
    void PrintOutResult();
    
    
    // in which interval is point r, return i for interval r_i r_{i+1}, -1 for out of range
    int getInterval(double &r);
    
protected:    
    // the grid points
    ub::vector<double> _r;
    
    // unknowns f_i and f"_i in one vector of 2*(n+1) size
    // they are stored in one array since they are used together 
    // in the linear equation system for the force matching scheme
    ub::vector<double> _f; 
    
    eBoundary _boundaries;
};

inline int CubicSpline::GenerateGrid(double min, double max, double h)
{
    int vec_size = 1 + (int)((max-min)/h);  //check it!
    _r.resize(vec_size);
    int i;
   
    double r_init;
   
    for (r_init = min, i=0; i < vec_size; r_init += h ) {
            _r[i++]= r_init;
    }
    _f.resize(2 * _r.size(), false);
    return _r.size();
}


inline double CubicSpline::A(double &r)
{
    return ( 1.0 - (r -_r[getInterval(r)])/(_r[getInterval(r)+1]-_r[getInterval(r)]) );
}

inline double CubicSpline::B(double &r)
{
    return  (r -_r[getInterval(r)])/(_r[getInterval(r)+1]-_r[getInterval(r)]) ;
}

inline double CubicSpline::C(double &r)
{
    double xxi, h;
    xxi = r -_r[getInterval(r)];
    h   = _r[getInterval(r)+1]-_r[getInterval(r)];
    
    return ( 0.5*xxi*xxi - (1.0/6.0)*xxi*xxi*xxi/h - (1.0/3.0)*xxi*h) ;
}

inline double CubicSpline::D(double &r)
{
    double xxi, h;
    xxi = r -_r[getInterval(r)];
    h   = _r[getInterval(r)+1]-_r[getInterval(r)]; 
    
    return ( (1.0/6.0)*xxi*xxi*xxi/h - (1.0/6.0)*xxi*h ) ;
}

inline int CubicSpline::getInterval(double &r)
{
    if (r < _r[0] || r > _r[_r.size() - 1]) return -1;
    return int( (r - _r[0]) / (_r[_r.size()-1] - _r[0]) * (_r.size() - 1) );
}

inline double CubicSpline::A_prime(int i, int flag)
{
    if (flag == 0) return -1.0/(_r[i+1]-_r[i]);
    else if (flag == 1) return -1.0/(_r[i+2]-_r[i+1]);
    else
        throw std::invalid_argument("Error: wrong flag in CubicSpline::A_prime"); 
}

inline double CubicSpline::B_prime(int i, int flag)
{
    if (flag == 0) return 1.0/(_r[i+1]-_r[i]);
    else if (flag == 1) return 1.0/(_r[i+2]-_r[i+1]);
    else
        throw std::invalid_argument("Error: wrong flag in CubicSpline::B_prime");
}

inline double CubicSpline::C_prime(int i, int flag)
{
    if (flag == 0) return (1.0/6.0)*(_r[i+1]-_r[i]);
    else if (flag == 1) return -(1.0/3.0)*(_r[i+2]-_r[i+1]);
    else
        throw std::invalid_argument("Error: wrong flag in CubicSpline::C_prime");  
}

inline double CubicSpline::D_prime(int i, int flag)
{
    if (flag == 0) return (1.0/3.0)*(_r[i+1]-_r[i]);
    else if (flag == 1) return -(1.0/6.0)*(_r[i+2]-_r[i+1]);
    else
        throw std::invalid_argument("Error: wrong flag in CubicSpline::D_prime");        
}

inline double CubicSpline::A_prime_l(int i)
{
    return -1.0/(_r[i+1]-_r[i]);
}

inline double CubicSpline::B_prime_l(int i)
{
    return 1.0/(_r[i+1]-_r[i]);
}

inline double CubicSpline::C_prime_l(int i)
{
    return (1.0/6.0)*(_r[i+1]-_r[i]);
}

inline double CubicSpline::D_prime_l(int i)
{
    return (1.0/3.0)*(_r[i+1]-_r[i]);
}

inline double CubicSpline::A_prime_r(int i)
{
    return -1.0/(_r[i+2]-_r[i+1]);
}

inline double CubicSpline::B_prime_r(int i)
{
    return 1.0/(_r[i+2]-_r[i+1]);
}

inline double CubicSpline::C_prime_r(int i)
{
    return -(1.0/3.0)*(_r[i+2]-_r[i+1]);
}

inline double CubicSpline::D_prime_r(int i)
{
    return -(1.0/6.0)*(_r[i+2]-_r[i+1]);
}

inline void CubicSpline::GetResult(ub::vector<double> *x_pointer) 
{
    ub::vector<double>::iterator ia;
    int i = 0;
    for (ia = x_pointer->begin(); ia != x_pointer->end(); ++ia, ++i) {
        _f[i] = *ia;
    }
}

inline double CubicSpline::getFunctionValue(double &r)
{
    int n = _f.size()/2;
    return -A(r)*_f[getInterval(r)] - B(r)*_f[getInterval(r) + 1] - C(r)*_f[n + getInterval(r)] - D(r)*_f[n + getInterval(r) + 1];
}

inline void CubicSpline::PrintOutResult()
{
    for (double x = _r[0]; x <= _r[_r.size() - 1]; x += 0.0001) {
        std::cout << x << " " << getFunctionValue(x) << "\n";
    }
}

#endif	/* _CUBICSPLINE_H */

