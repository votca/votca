
#include "potentialfunctioncbspl.h"

PotentialFunctionCBSPL::PotentialFunctionCBSPL(const int nlam_, const int ncutcoeff_,
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

    _dr = (_cut_off - _min )/( double (_nbreak - 1) );

    _M.resize(4,4,false);
    _M.clear();

    _M(0,0) =  1.0; _M(0,1) =  4.0; _M(0,2) =  1.0; _M(0,3) = 0.0;
    _M(1,0) = -3.0; _M(1,1) =  0.0; _M(1,2) =  3.0; _M(1,3) = 0.0;
    _M(2,0) =  3.0; _M(2,1) = -6.0; _M(2,2) =  3.0; _M(2,3) = 0.0;
    _M(3,0) = -1.0; _M(3,1) =  3.0; _M(3,2) = -3.0; _M(3,3) = 1.0;

    _M /= 6.0;

}

double PotentialFunctionCBSPL::CalculateF (const double r) const {

    if ( r < _min || r >= _cut_off ){

        return 0.0;

    }else {

        ub::vector<double> R;
        ub::vector<double> B;

        int indx = min( int( ( r - _min )/_dr ), _nbreak-2 );

        double rk = _min + indx*_dr;

        double t = ( r - rk)/_dr;

        R.resize(4,false); R.clear();

        R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;

        ub::vector<double> RM = ub::prod(R,_M);

        B.resize(4,false); B.clear();

        B(0) = _lam(indx); B(1) = _lam(indx+1); B(2) = _lam(indx+2);
        B(3) = _lam(indx+3);

        double u = ub::inner_prod(B,RM);

        return u;
    }

}
    // calculate first derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateDF(const int i, const double r) const{

    if ( r >= _min && r < _cut_off ) {

        int indx = min( int( ( r - _min )/_dr ), _nbreak-2 );

        if ( i >= indx && i <= indx+3 ){

            ub::vector<double> R;

            double rk = _min + indx*_dr;

            double t = ( r - rk)/_dr;

            R.resize(4,false); R.clear();

            R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;

            ub::vector<double> RM = ub::prod(R,_M);

            return RM(i-indx);

        }else{

            return 0.0;

        }


    } else {

        return 0;

    }
}
    // calculate second derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateD2F(const int i, const int j, const double r) const {

    // for cubic B-SPlines D2F is zero for all lamdas
    return 0.0;

}
