/* 
 * File:   CubicSpline.h
 * Author: ruehle
 *
 * Created on July 4, 2008, 3:03 PM
 */

#ifndef _CUBICSPLINE_H
#define	_CUBICSPLINE_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
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
    
    void setBondaries(eBoundary bc) {_boundaries = bc; }
    
    // Generates the r_k, returns the number of grid points
    int GenerateGrid(double min, double max, double h);
    // in which interval is point r, return i for interval r_i r_{i+1}, -1 for out of range
    int getInterval(double &r);

    // give string in rangeparser format: e.g. 1:0.1:10;11:1:20
    //int GenerateGrid(string range);

    // construct an interpolation spline
    void Interpolate(ub::vector<double> &x, ub::vector<double> &y);    
    // fit spline through data
    void Fit(ub::vector<double> &x, ub::vector<double> &y);
    
    
    double Calculate(double &x);
    
    template<typename vector_type1, typename vector_type2>
    void Calculate(vector_type1 &x, vector_type2 &y);
       
    // store data in the spline, todo: rename that function
    template<typename vector_type>
    void setSplineData(vector_type &f) { _f = f; }
    
    void Print(std::ostream &out, double interval = 0.0001 );
        
    ub::vector<double> &getX() {return _r; }
    ub::vector<double> &getSplineData() { return _f; }
    
    // stuff to construct fitting matrices
    
    // add the points
    template<typename matrix_type>
    void AddToFitMatrix(matrix_type &A, double x,
            int offset1, int offset2=0, double scale=1);
    
    template<typename matrix_type, typename vector_type>
    void AddToFitMatrix(matrix_type &M, vector_type &x, 
            int offset1, int offset2=0);
    
    // add the boundary conditions
    template<typename matrix_type>
    void AddBCToFitMatrix(matrix_type &A,
            int offset1, int offset2=0);


protected:    
    // the grid points
    ub::vector<double> _r;
    
    // unknowns f_i and f"_i in one vector of 2*(n+1) size
    // they are stored in one array since they are used together 
    // in the linear equation system for the force matching scheme
    ub::vector<double> _f; 
    
    eBoundary _boundaries;
    
    // A spline can be written in the form
    // S_i(x) =   A(x,x_i,x_i+1)*f_i     + B(x,x_i,x_i+1)*f''_i 
    //          + C(x,x_i,x_i+1)*f_{i+1} + D(x,x_i,x_i+1)*f''_{i+1}
    double A(double &r);
    double B(double &r);
    double C(double &r);
    double D(double &r);
  
    // tabulated derivatives at grid points. Second argument: 0 - left, 1 - right    
    double A_prime_l(int i);     
    double A_prime_r(int i);     
    double B_prime_l(int i);    
    double B_prime_r(int i);    
    double C_prime_l(int i);
    double C_prime_r(int i);
    double D_prime_l(int i);
    double D_prime_r(int i);
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

inline double CubicSpline::Calculate(double &r)
{
    int n = _f.size()/2;
    return -A(r)*_f[getInterval(r)] - B(r)*_f[getInterval(r) + 1] - C(r)*_f[n + getInterval(r)] - D(r)*_f[n + getInterval(r) + 1];
}

template<typename vector_type1, typename vector_type2>
inline void CubicSpline::Calculate(vector_type1 &x, vector_type2 &y)
{
    int n = _r.size();
    for(int i=0; i<x.size(); ++i)
        y(i) =  - A(x(i))*_f[getInterval(x(i))] 
                - B(x(i))*_f[getInterval(x(i)) + 1] 
                - C(x(i))*_f[n + getInterval(x(i))] 
                - D(x(i))*_f[n + getInterval(x(i)) + 1];
}

inline void CubicSpline::Print(std::ostream &out, double interval)
{
    for (double x = _r[0]; x <= _r[_r.size() - 1]; x += interval)
        out << x << " " << Calculate(x) << "\n";    
}

template<typename matrix_type>
inline void CubicSpline::AddToFitMatrix(matrix_type &M, double x, 
            int offset1, int offset2, double scale)
{
    int spi = getInterval(x);
    M(offset1, offset2 + spi) = A(x)*scale;
    M(offset1, offset2 + spi+1) = B(x)*scale;
    M(offset1, offset2 + spi + _r.size()) = C(x)*scale;
    M(offset1, offset2 + spi + _r.size() + 1) = D(x)*scale;
}

template<typename matrix_type, typename vector_type>
inline void CubicSpline::AddToFitMatrix(matrix_type &M, vector_type &x, 
            int offset1, int offset2)
{
    for(int i=0; i<x.size(); ++i) {
        int spi = getInterval(x(i));
        M(offset1+i, offset2 + spi) = A(x(i));
        M(offset1+i, offset2 + spi+1) = B(x(i));
        M(offset1+i, offset2 + spi + _r.size()) = C(x(i));
        M(offset1+i, offset2 + spi + _r.size() + 1) = D(x(i));
    }
}

template<typename matrix_type>
inline void CubicSpline::AddBCToFitMatrix(matrix_type &M,
            int offset1, int offset2)
{
    for(int i=0; i<_r.size() - 2; ++i) {
            M(offset1+i+1, offset2 + i) = A_prime_l(i);
            M(offset1+i+1, offset2 + i+1) = B_prime_l(i) - A_prime_r(i);
            M(offset1+i+1, offset2 + i+2) = -B_prime_r(i);

            M(offset1+i+1, offset2 + _r.size() + i) = C_prime_l(i);
            M(offset1+i+1, offset2 + _r.size() + i+1) = D_prime_l(i) - C_prime_r(i);
            M(offset1+i+1, offset2 + _r.size() + i+2) = -D_prime_r(i);
    }
    // currently only natural boundary conditions:
    M(offset1, offset2 + _r.size()) = 1;
    M(offset1 + _r.size() - 1, offset2 + 2*_r.size()-1) = 1;
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

#endif	/* _CUBICSPLINE_H */

