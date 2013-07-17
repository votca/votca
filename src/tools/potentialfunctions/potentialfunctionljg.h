/* 
 * 
 * Author: mashaya1
 *
 * Created on November 15, 2011, 9:51 PM
 */

#ifndef POTENTIALFUNCTIONLJG_H
#define	POTENTIALFUNCTIONLJG_H
#include "potentialfunction.h"

// LJ 12-6 potential class
// with c12,c6 parameters
class PotentialFunctionLJG : public PotentialFunction {
public:
    PotentialFunctionLJG(const string& name_,const double min_ = 0.0,
            const double max_ = 10.0);
    ~PotentialFunctionLJG() {};
    // calculate function value for given r
    double CalculateF (const double r) const;
    // calculate first derivative w.r.t. ith parameter
    double CalculateDF(const int i, const double r) const;
    // calculate second derivative w.r.t. ith parameter
    double CalculateD2F(const int i, const int j, const double r) const;
};




#endif	/* POTFUNCTION_LJG_H */

