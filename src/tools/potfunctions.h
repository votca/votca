/* 
 * File:   potfunctions.h
 * Author: sikandar
 *
 * Created on November 8, 2011, 11:52 PM
 */

#ifndef POTFUNCTIONS_H
#define	POTFUNCTIONS_H

#include <votca/tools/table.h>
#include <boost/numeric/ublas/vector.hpp>
#include <math.h>
using namespace std;
//namespace ub = boost::numeric::ublas;
using namespace votca::tools;

class FunctionForm {
public:
    FunctionForm() {}

    virtual ~FunctionForm() {}
    // read parameters from the input file
    void setParam(string filename) {

        Table param;
        param.Load(filename);

        if( param.size() != _nlam) {
            // need to throw better error
            throw std::runtime_error("parameters size mismatch");
        } else {
            for( int i = 0; i < _nlam; i++)
                _lam(i) = param.y(i);
        }
    }
    // set all parameters
    void setParam(const ub::vector<double> param) { _lam = param; _nlam = param.size(); }
    // set ith parameter
    void setParam(const int i, const double val) { _lam[i] = val; }
    // set minimum r value to avoid large values
    void setMinDist(const double min) { _min = min; }
    // set cut-off value
    void setCutOffDist(const double cutoff) { _cut_off = cutoff; }
    // calculate function
    virtual double CalculateF (const double r) const = 0;
    // calculate first derivative w.r.t. ith parameter
    virtual double CalculateDF(const int i, const double r) const = 0;
    // calculate second derivative w.r.t. ith parameter
    virtual double CalculateD2F(const int i, const int j, const double r) const = 0;
    // calculate integrant r.f(r)
    virtual double CalculateIntRF(const double r) const = 0;
    // return parameter
    ub::vector<double> getParam() const { return _lam; }
    // return ith parameter
    double getParam(const int i) const { return _lam[i]; }
    // return size of parameters
    int getParamSize() const { return _nlam; }
    // return cut-off value
    double getCutOff() const { return _cut_off; }

protected:
    ub::vector<double> _lam;
    int _nlam;
    double _cut_off;
    double _min;
};

// LJ 12-6 potential class

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
 
double FunctionLJ126::CalculateF (const double r) const {

    if ( r >= _min && r <= _cut_off ) {

        return 4.0*_lam(1)*( pow((_lam(0)/r),12) - pow((_lam(0)/r),6) ) ;

    } else {

        return 0.0;

    }
}
    // calculate first derivative w.r.t. ith parameter
double FunctionLJ126::CalculateDF(const int i, const double r) const {

    if ( r >= _min && r <= _cut_off ) {

        switch(i) {
            case 0:
                return 4.0*_lam(1)*( 12.0*pow(_lam(0),11)/pow(r,12) - 6.0*pow(_lam(0),5)/pow(r,6) );
            case 1:
                return 4.0*( pow((_lam(0)/r),12) - pow((_lam(0)/r),6) );
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
                        return 4.0*_lam(1)*( 132.0*pow(_lam(0),10)/pow(r,12)
                                            - 30.0*pow(_lam(0),4)/pow(r,6) );
                    case 1:
                        return 4.0*( 12.0*pow(_lam(0),11)/pow(r,12)
                                     - 6.0*pow(_lam(0),5)/pow(r,6) );
                }
                
            case 1:
                switch(j){
                    case 0:
                        return 4.0*( 12.0*pow(_lam(0),11)/pow(r,12)
                                     - 6.0*pow(_lam(0),5)/pow(r,6) );
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
    return 4.0*_lam(1)*( pow(_lam(0),12)/(-11.0*pow(r,11)) - pow(_lam(0),6)/(-5.0*pow(r,5)) );
}

#endif