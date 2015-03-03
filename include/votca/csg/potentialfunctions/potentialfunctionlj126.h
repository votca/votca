/*
 * Author: sikandar
 *
 * Created on November 8, 2011, 11:52 PM
 */

#ifndef POTENTIALFUNCTIONLJ126_H
#define	POTENTIALFUNCTIONLJ126_H

#include <votca/tools/table.h>
#include <boost/numeric/ublas/vector.hpp>
#include <math.h>
#include "potentialfunction.h"

using namespace std;
using namespace votca::tools;

// LJ 12-6 potential class
// with c12,c6 parameters
class PotentialFunctionLJ126 : public PotentialFunction {
public:
    PotentialFunctionLJ126(const string& name_,const double min_=0.0,
	const double max_=10.0);
    ~PotentialFunctionLJ126(){};
    // calculate function value for given r
    double CalculateF (const double r) const;
    // calculate first derivative w.r.t. ith parameter
    double CalculateDF(const int i, const double r) const;
    // calculate second derivative w.r.t. ith parameter
    double CalculateD2F(const int i, const int j, const double r) const;

};

#endif
