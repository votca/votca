
#include "potentialfunction.h"

void PotentialFunction::setParam(string filename) {

        Table param;
        param.Load(filename);

        if( param.size() != _nlam) {

            throw std::runtime_error("Potential parameters size mismatch!\n"
                    "Check input parameter file \""
                    + filename + "\" \nThere should be "
                    + boost::lexical_cast<string>( _nlam ) + " parameters");
        } else {
            for( int i = 0; i < _nlam; i++)
                _lam(i) = param.y(i);
            
        }

    }
