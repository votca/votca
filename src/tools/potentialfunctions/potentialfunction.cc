
#include "potentialfunction.h"

PotentialFunction::PotentialFunction(const int nlam_,
		        const double min_,const double max_){
   _lam.resize(nlam_);
   _lam.clear();
   _min = min_;
   _cut_off = max_;	

}


void PotentialFunction::setParam(string filename) {

        Table param;
        param.Load(filename);

        if( param.size() != _lam.size()) {

            throw std::runtime_error("Potential parameters size mismatch!\n"
                    "Check input parameter file \""
                    + filename + "\" \nThere should be "
                    + boost::lexical_cast<string>( _lam.size() ) + " parameters");
        } else {
            for( int i = 0; i < _lam.size(); i++)
                _lam(i) = param.y(i);
            
        }

    }

void PotentialFunction::SaveParam(const string& filename){

    Table param;
    param.SetHasYErr(false);
    param.resize(_lam.size(), false);

    for (int i = 0; i < _lam.size(); i++){
        param.set(i, i, _lam(i), 'i');
    }
    param.Save(filename);

}
