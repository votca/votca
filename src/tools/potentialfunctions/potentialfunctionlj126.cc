
#include "potentialfunctionlj126.h"

PotentialFunctionLJ126::PotentialFunctionLJ126() {

    _nlam = 2;
    _lam.resize(_nlam);
    _lam.clear();
    _min = 0.0; // default rmin = 0.0
    _cut_off = 100.0; // default cut-off 100 nm.

}

PotentialFunctionLJ126::PotentialFunctionLJ126(const double min_, const double max_){

   _nlam = 2;
   _lam.resize(_nlam);
   _lam.clear();
   _min = min_;
   _cut_off = max_;

}

double PotentialFunctionLJ126::CalculateF (const double r) const {

    if ( r >= _min && r <= _cut_off ) {

        return _lam(0)/pow(r,12) - _lam(1)/pow(r,6) ;

    } else {

        return 0.0;

    }
}
    // calculate first derivative w.r.t. ith parameter
double PotentialFunctionLJ126::CalculateDF(const int i, const double r) const {

    if ( r >= _min && r <= _cut_off ) {

        switch(i) {
            case 0:
                return 1.0/pow(r,12);
            case 1:
                return -1.0/pow(r,6);
        }

    } else {

        return 0;

    }
}
    // calculate second derivative w.r.t. ith parameter
double PotentialFunctionLJ126::CalculateD2F(const int i, const int j, const double r) const {

    if ( r >= _min && r <= _cut_off ) {
        switch(i) {
            case 0:
                switch(j){
                    case 0:
                        return 0.0;
                    case 1:
                        return 0.0;
                }

            case 1:
                switch(j){
                    case 0:
                        return 0.0;
                    case 1:
                        return 0.0;
                }
        }

    } else {

        return 0;

    }
}
