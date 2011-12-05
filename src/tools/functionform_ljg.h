/* 
 * File:   potfunction_LJG.h
 * Author: mashaya1
 *
 * Created on November 15, 2011, 9:51 PM
 */

#ifndef POTFUNCTION_LJG_H
#define	POTFUNCTION_LJG_H
#include "functionform.h"

// LJ 12-6 potential class
// with c12,c6 parameters
class FunctionLJG : public FunctionForm {
public:
    FunctionLJG();
    FunctionLJG(const double min_, const double max_);
    ~FunctionLJG() {};
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

FunctionLJG::FunctionLJG() {

    _nlam = 5;
    _lam.resize(_nlam);
    _lam.clear();
    _min = 0.0;
    _cut_off = 100; // default cut-off 100 nm.
}

FunctionLJG::FunctionLJG(const double min_, const double max_){

   _nlam = 5;
   _lam.resize(_nlam);
   _lam.clear();
   _min = min_;
   _cut_off = max_;

}

double FunctionLJG::CalculateF (const double r) const {


    if ( r >= _min && r <= _cut_off ) {

        // lj 12-6 part + gaussian
        return _lam(0)/pow(r,12) - _lam(1)/pow(r,6)
                + _lam(2)*exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );

    } else {

        return 0.0;

    }

}
    // calculate first derivative w.r.t. ith parameter
double FunctionLJG::CalculateDF(const int i, const double r) const {
    
    if ( r >= _min && r <= _cut_off ) {

        switch(i) {
            case 0:
                return 1.0/pow(r,12);
            case 1:
                return -1.0/pow(r,6);
            case 2:
                return exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
            case 3:
                return -1.0*_lam(2)*(r-_lam(4))*(r-_lam(4)) *
                       exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
            case 4:
                return 2.0*_lam(2)*_lam(3)*(r-_lam(4)) *
                       exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );                
        }

    } else {

        return 0.0;

    }
    
}
    // calculate second derivative w.r.t. ith parameter
double FunctionLJG::CalculateD2F(const int i, const int j, const double r) const {

    if ( r >= _min && r <= _cut_off ) {
        
        switch(i) {
            case 0:
                // all second derivatives w.r.t. c12 are zero
                return 0.0;
            case 1:
                return 0.0;
            case 2:
                switch(j){
                    case 0:
                        return 0.0;
                    case 1:
                        return 0.0;
                    case 2:
                        return 0.0;
                    case 3:
                        return -1.0*(r-_lam(4))*(r-_lam(4)) *
                                exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
                    case 4:
                        return 2.0*_lam(3)*(r-_lam(4)) *
                                exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
                }                                
            case 3:
                switch(j){
                    case 0:
                        return 0.0;
                    case 1:
                        return 0.0;
                    case 2:
                        return -1.0*(r-_lam(4))*(r-_lam(4)) *
                                exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
                    case 3:
                        return _lam(2)*pow((r-_lam(4)),4) *
                                exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
                    case 4:
                        return 2.0*_lam(2)*(r-_lam(4)) *
                                ( 1.0 - _lam(3)*pow((r-_lam(4)),2) ) *
                                exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
                } 
            case 4:
                switch(j){
                    case 0:
                        return 0.0;
                    case 1:
                        return 0.0;
                    case 2:
                        return 2.0*_lam(3)*(r-_lam(4)) *
                                exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
                    case 3:
                        return 2.0*_lam(2)*(r-_lam(4)) *
                                ( 1.0 - _lam(3)*pow((r-_lam(4)),2) ) *
                                exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
                        
                    case 4:
                        return 2.0*_lam(2)*_lam(3)*
                                ( 2.0*_lam(3)*pow((r-_lam(4)),2) - 1.0 ) *
                                exp( -1.0*_lam(3)*(r-_lam(4))*(r-_lam(4)) );
                } 
        }
    } else {

        return 0.0;

    }    
}
    // calculate integrant r.f(r)
double FunctionLJG::CalculateIntRF(const double r) const{
    return 0.0;
}


#endif	/* POTFUNCTION_LJG_H */

