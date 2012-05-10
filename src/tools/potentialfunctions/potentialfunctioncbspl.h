
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
    PotentialFunctionCBSPL(const string& name_,const int nlam_,const int ncutcoeff_,
     const double min_=0.0, const double max_=10.0);
    ~PotentialFunctionCBSPL(){}
    // calculate function value for given r
    double CalculateF (const double r) const;
    // calculate first derivative w.r.t. ith parameter
    double CalculateDF(const int i, const double r) const;
    // calculate second derivative w.r.t. ith parameter
    double CalculateD2F(const int i, const int j, const double r) const;
    
    int getOptParamSize() const ;
    
    void setParam(string filename);
    // save parameters to the file
    void SaveParam(const string& filename);
    
    void SavePotTab(const string& filename, const double step);

    void SavePotTab(const string& filename, const double step, const double rmin, const double rcut);    
    void setOptParam(const int i, const double val);
    
    double getOptParam(const int i) const;
    
    void extrapolExclParam();

protected:

    // exclude these many first coefficients from optimization
    // since the region relevant to these coefficients is not sampled
    // the value of _nexcl is determined from rmin
    int _nexcl;
    // fix these many coeff near the cut-off to zero to ensure
    // zero potential and force values near cut-off
    int _ncutcoeff;

    int _nbreak;
    double _dr;
    ub::vector<double> _rbreak;
    
    ub::matrix<double> _M;
    
    

};

#endif
