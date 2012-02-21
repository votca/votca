
#ifndef POTENTIALFUNCTIONCBSPL_H
#define	POTENTIALFUNCTIONCBSPL_H

#include <votca/tools/table.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <math.h>
#include "potentialfunction.h"

using namespace std;
using namespace votca::tools;

class PotentialFunctionCBSPL : public PotentialFunction {
public:
    PotentialFunctionCBSPL(const int nlam_, const int ncutcoeff_,
                           const double min_, const double max_);
    ~PotentialFunctionCBSPL(){}
    // calculate function value for given r
    double CalculateF (const double r) const;
    // calculate first derivative w.r.t. ith parameter
    double CalculateDF(const int i, const double r) const;
    // calculate second derivative w.r.t. ith parameter
    double CalculateD2F(const int i, const int j, const double r) const;

protected:

    // fix these many coeff near the cut-off to zero to ensure
    // zero potential and force values near cut-off
    int _ncutcoeff;
    // total coeff = _nlam + _nceoff
    int _ncoeff;
    // total break points
    int _nbreak;

    double _dr;

    ub::matrix<double> _M;

};

#endif
