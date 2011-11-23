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
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>

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
            
            throw std::runtime_error("Potential parameters size mismatch!\nCheck input parameter file \""
                    + filename + "\" \nThere should be " + boost::lexical_cast<string>( _nlam ) + " parameters");
        } else {
            for( int i = 0; i < _nlam; i++)
                _lam(i) = param.y(i);
            
        }

    }
    // set all parameters
    void setParam(const ub::vector<double> param) { _lam = param; _nlam = param.size(); }
    // set ith parameter
    void setParam(const int i, const double val) { _lam(i) = val; }
    // fit parameters to given potential
    virtual void fitParam(string filename) = 0;
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
    double getParam(const int i) const { return _lam(i); }
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

#endif
