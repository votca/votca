/* 
 * File:   potfunctions.h
 * Author: sikandar
 *
 * Created on November 8, 2011, 11:52 PM
 */

#ifndef POTFUNCTION_LJ126_H
#define	POTFUNCTION_LJ126_H

#include <votca/tools/table.h>
#include <boost/numeric/ublas/vector.hpp>
#include <math.h>
#include "potfunctions.h"
using namespace std;
//namespace ub = boost::numeric::ublas;
using namespace votca::tools;

// LJ 12-6 potential class
// with c12,c6 parameters
class FunctionLJ126 : public FunctionForm {
public:
    FunctionLJ126();
    FunctionLJ126(const double min_, const double max_);
    ~FunctionLJ126(){};
    // fit parameters to given potential
    void fitParam(string filename) {}
    // calculate function value for given r
    double CalculateF (const double r) const;
    // calculate first derivative w.r.t. ith parameter
    double CalculateDF(const int i, const double r) const;
    // calculate second derivative w.r.t. ith parameter
    double CalculateD2F(const int i, const int j, const double r) const;
    // calculate integrant r.f(r)
    double CalculateIntRF(const double r) const;
};

FunctionLJ126::FunctionLJ126() {

    _nlam = 2;
    _lam.resize(_nlam);
    _lam.clear();
    _min = 0.0; // default rmin = 0.0
    _cut_off = 100.0; // default cut-off 100 nm.
    
}

FunctionLJ126::FunctionLJ126(const double min_, const double max_){

   _nlam = 2;
   _lam.resize(_nlam);
   _lam.clear();
   _min = min_;
   _cut_off = max_;
   
}

double FunctionLJ126::CalculateF (const double r) const {

    if ( r >= _min && r <= _cut_off ) {

        return _lam(0)/pow(r,12) - _lam(1)/pow(r,6) ;

    } else {

        return 0.0;

    }
}
    // calculate first derivative w.r.t. ith parameter
double FunctionLJ126::CalculateDF(const int i, const double r) const {

    if ( r >= _min && r <= _cut_off ) {

        switch(i) {
            case 0:
                return 1.0/pow(r,12);
            case 1:
                return -1.0/pow(r,6);
        }

    } else {

        return 0;

    }
}
    // calculate second derivative w.r.t. ith parameter
double FunctionLJ126::CalculateD2F(const int i, const int j, const double r) const {

    if ( r >= _min && r <= _cut_off ) {
        switch(i) {
            case 0:
                switch(j){
                    case 0:
                        return 0.0;
                    case 1:
                        return 0.0;
                }

            case 1:
                switch(j){
                    case 0:
                        return 0.0;
                    case 1:
                        return 0.0;
                }
        }

    } else {

        return 0;

    }
}
    // calculate integrant r.f(r)
double FunctionLJ126::CalculateIntRF(const double r) const{
    return 0.0;
}



#endif