/* 
 * File:   functionform.h
 * Author: sikandar
 *
 * Created on October 22, 2011, 7:33 PM
 */

#include<boost/numeric/ublas/vector.hpp>
#include <votca/tools/table.h>
#include<iostream>
//namespace ub = boost::numeric::ublas;
using namespace votca::tools;

class FunctionForm {
public:
    FunctionForm() {}
    
    virtual ~FunctionForm() {}
    // read parameters from the input file
    void setParam(string filename) {}
    // set all parameters
    void setParam(const ub::vector<double> param) { _lam = param; _nlam = param.size(); }
    // set ith parameter
    void setParam(const int i, const double val) { _lam[i] = val; }
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
    int getParamSize() { return _nlam; }
    
    
protected:
    ub::vector<double> _lam;
    int _nlam;
};
