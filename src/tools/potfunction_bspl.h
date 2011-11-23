/* 
 * File:   potfunctions.h
 * Author: sikandar
 *
 * Created on November 8, 2011, 11:52 PM
 */

#ifndef POTFUNCTION_BSPL_H
#define	POTFUNCTION_BSPL_H

#include <votca/tools/table.h>
#include <boost/numeric/ublas/vector.hpp>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include "potfunctions.h"

using namespace std;
//namespace ub = boost::numeric::ublas;
using namespace votca::tools;

// LJ 12-6 potential class
// with c12,c6 parameters
class FunctionBSPL : public FunctionForm {
public:
    FunctionBSPL(const int nlam_, const int ncutcoeff_,
                           const double min_, const double max_);
    ~FunctionBSPL();
    // fit parameters to given potential
    void fitParam(string filename);
    // calculate function value for given r
    double CalculateF (const double r) const;
    // calculate first derivative w.r.t. ith parameter
    double CalculateDF(const int i, const double r) const;
    // calculate second derivative w.r.t. ith parameter
    double CalculateD2F(const int i, const int j, const double r) const;
    // calculate integrant r.f(r)
    double CalculateIntRF(const double r) const;

protected:

    // fix these many coeff near the cut-off to zero to ensure
    // zero potential and force values near cut-off
    int _ncutcoeff;
    // total coeff = _nlam + _nceoff
    int _ncoeff;
    // total break points
    int _nbreak;

    // GSL stuff for B-splines
    gsl_bspline_workspace *_gsl_bw;
    
    
};

// 
FunctionBSPL::FunctionBSPL(const int nlam_, const int ncutcoeff_,
                           const double min_, const double max_) {

    /* Here _nlam is the total number of coeff values that are to be optimized
     * To ensure that potential and force go to zero smoothly near cut-off,
     * as suggested in Ref. PCCP, 11, 1901, 2009, coeff values leading up to
     * cut-off and beyond take a value of zero.
     * Hence, total coeff = _nlam + _ncutcoeff ( i.e. _nlam variable coeff
     *  and _ncutcoeff zero coeff.
     */

    _nlam = nlam_;

    _ncutcoeff = ncutcoeff_;

    _ncoeff = _nlam + _ncutcoeff;

    _lam.resize(_ncoeff,false);
    
    _lam.clear();

    _nbreak = _ncoeff - 2;

    _min = min_;
    _cut_off = max_;

    /* allocate a cubic bspline workspace (k = 4) */
    _gsl_bw = gsl_bspline_alloc(4, _nbreak);
    /* use uniform breakpoints on [_min,_cut_off] */
    gsl_bspline_knots_uniform(_min, _cut_off, _gsl_bw);
    
}

double FunctionBSPL::CalculateF (const double r) const {

    if ( r >= _min && r < _cut_off ) {        

        gsl_vector *B;
        B = gsl_vector_alloc(_ncoeff);
        gsl_bspline_eval(r, B, _gsl_bw);

        double u = 0.0;
        for( int row = 0; row < _ncoeff; row++)
            u += _lam(row) * gsl_vector_get(B, row);

        gsl_vector_free(B);
        return u;

    } else {

        return 0.0;

    }
}
    // calculate first derivative w.r.t. ith parameter
double FunctionBSPL::CalculateDF(const int i, const double r) const{

    if ( r >= _min && r < _cut_off ) {

        gsl_vector *B;
        B = gsl_vector_alloc(_ncoeff);
        gsl_bspline_eval(r, B, _gsl_bw);

        double u = gsl_vector_get(B,i);

        gsl_vector_free(B);

        return u;


    } else {

        return 0;

    }
}
    // calculate second derivative w.r.t. ith parameter
double FunctionBSPL::CalculateD2F(const int i, const int j, const double r) const {

    // for cubic B-SPlines D2F is zero for all lamdas
    return 0.0;

}
    // calculate integrant r.f(r)
double FunctionBSPL::CalculateIntRF(const double r) const {
    return 0.0;
}

void FunctionBSPL::fitParam(string filename){

    Table pot;
    pot.Load(filename);

    int n = pot.size();

    gsl_vector *y,*c;
    gsl_matrix *X, *cov;
    gsl_multifit_linear_workspace *mw;
    double chisq;
    gsl_vector *B;

    B = gsl_vector_alloc(_ncoeff);
    y = gsl_vector_alloc(n);
    c = gsl_vector_alloc(_ncoeff);
    X = gsl_matrix_alloc(n, _ncoeff);    
    cov = gsl_matrix_alloc(_ncoeff, _ncoeff);
    mw = gsl_multifit_linear_alloc(n, _ncoeff);
    
    /* construct the fit matrix X */
    for (int i = 0; i < n; ++i){

        double xi = pot.x(i);

        /* compute all B-splines for xi */
        gsl_bspline_eval(xi, B, _gsl_bw);

        /* fill in row i of X */
        for (int j = 0; j < _ncoeff; ++j){

            double Bj = gsl_vector_get(B, j);
            gsl_matrix_set(X, i, j, Bj);

        }
        // also store pot values in gsl_vector y
        gsl_vector_set(y, i, pot.y(i));

    }

    /* do the fit */
    gsl_multifit_linear(X, y, c, cov, &chisq, mw);

    /* copy first _nlam coefficients from c to _lam since last _ncutcoeff are zero
     * _lam was initially zero so NO need to set last _ncutcoeff _lam again zero
     */
    for( int i = 0; i < _nlam; i++){

        _lam(i) = gsl_vector_get(c,i);
        
    }
    
    cout << "Cubic B-spline fit quality:" << endl;
    cout << "\t \t" << "chisq = "<< chisq << endl;

    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_vector_free(B);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(mw);

    return ;
}

FunctionBSPL::~FunctionBSPL(){

     gsl_bspline_free(_gsl_bw);
    
}

#endif
