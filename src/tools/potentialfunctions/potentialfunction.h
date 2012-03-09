/* 
 * 
 * Author: sikandar
 *
 * Created on November 8, 2011, 11:52 PM
 */

#ifndef POTENTIALFUNCTION_H
#define	POTENTIALFUNCTION_H

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
using namespace votca::tools;

class PotentialFunction {
public:
    
    virtual ~PotentialFunction() {}
    // read parameters from the input file
    virtual void setParam(string filename);
    // save parameters to the file
    virtual void SaveParam(const string& filename);
    // write potential table
    virtual void SavePotTab(const string& filename, const double step);
    // set all parameters
    void setParam(const ub::vector<double> param){ _lam = param; }
    // set ith parameter
    void setParam(const int i, const double val) { _lam(i) = val; }
    // set ith parameter among those to be optimized
    virtual void setOptParam(const int i, const double val) { 
        setParam(i,val);
    }
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
    // return parameter
    ub::vector<double>& Params() { return _lam; }
    // return ith parameter
    double getParam(const int i) const { return _lam(i); }
    // return ith parameter among those to be optimized
    virtual double getOptParam(const int i) const {
        return getParam(i);
    }
    // return size of parameters
    int getParamSize() const { return _lam.size(); }
    // return size of parameters to be optimized
    virtual int getOptParamSize() const { return getParamSize();}
    // return cut-off value
    double getCutOff() const { return _cut_off; }
    double getMinDist() const { return _min; }

protected:
    
    PotentialFunction(const int nlam_,const double min_,const double max_);
    
    ub::vector<double> _lam;
    double _cut_off;
    double _min;
};

#endif
