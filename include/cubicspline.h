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
    
    // Generates the r_k, returns the number of grid points
    int GenerateGrid(double min, double max, double h);
    
    double A(double &r);
    double B(double &r);
    double C(double &r);
    double D(double &r);
  
    // tabulated derivatives at grid points. Second argument: 0 - left, 1 - right
    double A_prime(int i, int flag); 
    double B_prime(int i, int flag);    
    double C_prime(int i, int flag);
    double D_prime(int i, int flag);
    
    void GetResult(ub::vector<double> *x_pointer);
    double getFunctionValue(double &r);
    void PrintOutResult();
    
    
    // in which interval is point r, return i for interval r_i r_{i+1}, -1 for out of range
    inline int getInterval(double &r);
    
protected:
    

    ub::vector<double> _r;
    ub::vector<double> _f; // unknowns f_i and f"_i in one vector of 2*(n+1) size
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
    else {
        cout << "Error: wrong flag in CubicSpline::A_prime";
    }
    
}

inline double CubicSpline::B_prime(int i, int flag)
{
    if (flag == 0) return 1.0/(_r[i+1]-_r[i]);
    else if (flag == 1) return 1.0/(_r[i+2]-_r[i+1]);
    else {
        cout << "Error: wrong flag in CubicSpline::B_prime";
    }    
}

inline double CubicSpline::C_prime(int i, int flag)
{
    if (flag == 0) return (1.0/6.0)*(_r[i+1]-_r[i]);
    else if (flag == 1) return -(1.0/3.0)*(_r[i+2]-_r[i+1]);
    else {
        cout << "Error: wrong flag in CubicSpline::C_prime";
    }    
}

inline double CubicSpline::D_prime(int i, int flag)
{
    if (flag == 0) return (1.0/3.0)*(_r[i+1]-_r[i]);
    else if (flag == 1) return -(1.0/6.0)*(_r[i+2]-_r[i+1]);
    else {
        cout << "Error: wrong flag in CubicSpline::D_prime";
    }    
}

void CubicSpline::GetResult(ub::vector<double> *x_pointer) 
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
    return A(r)*_f[getInterval(r)] + B(r)*_f[getInterval(r) + 1] + C(r)*_f[n + getInterval(r)] + D(r)*_f[n + getInterval(r) + 1];
}

void CubicSpline::PrintOutResult()
{
    for (double x = _r[0]; x <= _r[_r.size() - 1]; x += 0.005) {
        cout << x << " " << getFunctionValue(x) << "\n";
    }
}
#endif	/* _CUBICSPLINE_H */

