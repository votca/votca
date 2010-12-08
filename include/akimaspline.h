/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _AKIMASPLINE_H
#define	_AKIMASPLINE_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <iostream>

namespace votca { namespace tools {

namespace ub = boost::numeric::ublas;

/**
    \brief An Akima Spline Class
 *
 *  does Akima interpolation based on "A new method of interpolation and smooth
 *  curve fitting based on local procedures"
 *  fit not yet included
 */

class AkimaSpline
{    
public:
    // default constructor
    AkimaSpline() :
        _boundaries(splineNormal) {}
    
    // destructor
    ~AkimaSpline() {};
    
    /// enum for type of boundary condition
    enum eBoundary {
        splineNormal = 0,  ///< normal boundary conditions: \f$f_0=f_N=0\f$ 
        splinePeriodic    ///< periodic boundary conditions: \f$f_0=f_N\f$
    };
    
    /// set the boundary type of the spline
    void setBC(eBoundary bc) {_boundaries = bc; }
    
    /// Generates the r_k, returns the number of grid points
    /// max is included in the interval.
    int GenerateGrid(double min, double max, double h);
    
    /// determine the interval the point r is in
    /// returns i for interval r_i r_{i+1}, -1 for out of range
    int getInterval(const double &r);

    /// getSlope: handles all special cases to determine the slope t based
    /// on m1,m2,m3,m4
    double getSlope(double m1, double m2, double m3, double m4);
    
    /// ERROR-PRONE implementation, make it better!!!
    double getGridPoint(const int &i) {return _r[i];}

    // give string in rangeparser format: e.g. 1:0.1:10;11:1:20
    //int GenerateGrid(string range);

    /// \brief construct an interpolation spline
    ///
    ///   x, y are the the points to construct interpolation,
    /// both vectors must be of same size
    void Interpolate(ub::vector<double> &x, ub::vector<double> &y);
    
    /// \brief fit spline through noisy data
    ///
    /// x,y are arrays with noisy data, both vectors must be of same size
    void Fit(ub::vector<double> &x, ub::vector<double> &y);
    
    
    /// Calculate the function value
    double Calculate(const double &x);

    /// Calculate the function derivative
    double CalculateDerivative(const double &x);
    
    
    /// Calculate the function value for a whole array, story it in y
    template<typename vector_type1, typename vector_type2>
    void Calculate(vector_type1 &x, vector_type2 &y);
    template<typename vector_type1, typename vector_type2>
    void CalculateDerivative(vector_type1 &x, vector_type2 &y);
       
    /// print out results
    void Print(std::ostream &out, double interval = 0.0001 );
        
    /// get the grid array x
    ub::vector<double> &getX() {return _r; }
    
protected:    
    // the grid points
    ub::vector<double> _r;
    
    eBoundary _boundaries;

    // p1,p2,p3,p4 and t1,t2 in Akima paper, page 591
    ub::vector<double> p0;
    ub::vector<double> p1;
    ub::vector<double> p2;
    ub::vector<double> p3;
    ub::vector<double> t;
};

inline int AkimaSpline::GenerateGrid(double min, double max, double h)
{
    int vec_size = (int)((max-min)/h+1.00000001);
    _r.resize(vec_size);
    int i;
   
    double r_init;
   
    for (r_init = min, i=0; i < vec_size - 1; r_init += h ) {
            _r[i++]= r_init;
    }
    _r[i] = max;
    return _r.size();
}

inline double AkimaSpline::Calculate(const double &r)
{
    int interval =  getInterval(r);
    double z = r-_r[interval];
    return  p0(interval)
            + p1(interval)*z
            + p2(interval)*z*z
            + p3(interval)*z*z*z;
}

inline double AkimaSpline::CalculateDerivative(const double &r)
{
    int interval =  getInterval(r);
    double z = r-_r[interval];
    return  + p1(interval)
            + 2.0*p2(interval)*z
            + 3.0*p3(interval)*z*z;
}

template<typename vector_type1, typename vector_type2>
inline void AkimaSpline::Calculate(vector_type1 &x, vector_type2 &y)
{
    for(size_t i=0; i<x.size(); ++i) 
        y(i) = Calculate(x(i));
}

template<typename vector_type1, typename vector_type2>
inline void AkimaSpline::CalculateDerivative(vector_type1 &x, vector_type2 &y)
{
    for(size_t i=0; i<x.size(); ++i)
        y(i) = CalculateDerivative(x(i));
}

inline void AkimaSpline::Print(std::ostream &out, double interval)
{
    for (double x = _r[0]; x < _r[_r.size() - 1]; x += interval)
        out << x << " " << Calculate(x) << "\n";    
}

inline int AkimaSpline::getInterval(const double &r)
{
    if (r < _r[0]) return 0;
    if(r > _r[_r.size() - 2]) return _r.size()-2;
    size_t i;
    for(i=0; i<_r.size(); ++i)
        if(_r[i]>r) break;
    return i-1;
}

inline double AkimaSpline::getSlope(double m1, double m2, double m3, double m4)
{
    if ((m1==m2) && (m3==m4)) {
        return (m2+m3)/2.0;
    } else {
        return (fabs(m4-m3)*m2 + fabs(m2-m1)*m3) / (fabs(m4-m3) + fabs(m2-m1));
    }
}

}}

#endif	/* _AKIMASPLINE_H */

