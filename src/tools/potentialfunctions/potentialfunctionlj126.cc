
#include "potentialfunctionlj126.h"

PotentialFunctionLJ126::PotentialFunctionLJ126(const string& name_,const double min_,
	const double max_) : PotentialFunction(name_,2,min_,max_){
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

        return 0.0;

    }
}

// calculate second derivative w.r.t. ith and jth parameters
double PotentialFunctionLJ126::CalculateD2F(const int i, const int j, 
        const double r) const {

    return 0.0;
    
}
