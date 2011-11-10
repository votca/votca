/* 
 * File:   functionlj126.h
 * Author: sikandar
 *
 * Created on October 23, 2011, 1:06 AM
 */

#ifndef FUNCTIONLJ126_H
#define	FUNCTIONLJ126_H

#include "functionform.h"
#include <boost/numeric/ublas/vector.hpp>
#include <math.h>
using namespace std;
//namespace ub = boost::numeric::ublas;
using namespace votca::tools;
class FunctionLJ126 : public FunctionForm {
public:
    FunctionLJ126() {_nlam = 2; _lam.resize(_nlam); _lam[0] = 0.317; _lam[1]=0.6503;}
    ~FunctionLJ126() {};
    // calculate function value for given r
    double CalculateF (const double r) const;
    // calculate first derivative w.r.t. ith parameter
    double CalculateDF(const int i, const double r) const;
    // calculate second derivative w.r.t. ith parameter
    double CalculateD2F(const int i, const int j, const double r) const;
    // calculate integrant r.f(r)
    double CalculateIntRF(const double r) const;
};

#endif	/* FUNCTIONLJ126_H */

double FunctionLJ126::CalculateF (const double r) const {
    return 4.0*_lam(1)*( pow((_lam(0)/r),12) - pow((_lam(0)/r),6) ) ;
}
    // calculate first derivative w.r.t. ith parameter
double FunctionLJ126::CalculateDF(const int i, const double r) const {
    switch(i) {
        case 1:
            return 4.0*_lam(1)*( 12.0*pow(_lam(0),11)/pow(r,12) - 6.0*pow(_lam(0),5)/pow(r,6) );
        case 2:
            return 4.0*( pow((_lam(0)/r),12) - pow((_lam(0)/r),6) );
    }
}
    // calculate second derivative w.r.t. ith parameter
double FunctionLJ126::CalculateD2F(const int i, const int j, const double r) const {
    switch(i) {
        case 1:
            return 4.0*_lam(1)*( 132.0*pow(_lam(0),10)/pow(r,12) - 30.0*pow(_lam(0),4)/pow(r,6) );
        case 2:
            return 0;
    }
}
    // calculate integrant r.f(r)
double FunctionLJ126::CalculateIntRF(const double r) const{
    return 4.0*_lam(1)*( pow(_lam(0),12)/(-11.0*pow(r,11)) - pow(_lam(0),6)/(-5.0*pow(r,5)) );
}