#ifndef FILE_CHARGES
#define FILE_CHARGES

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>

namespace votca { namespace moo {

using namespace std;

class multipoles
{
public: 
    vector<double> mpls; // bear charges for the molecules

    void cp_mpls( const multipoles &A)
    {
        if (mpls.size() != A.mpls.size() ){
            mpls.resize(A.mpls.size());
        }
        copy( A.mpls.begin(), A.mpls.end(), mpls.begin());
    }
	
    int read_crg_eps(const char *);

    const double & get_mpl(const int & i) const {
	return mpls[i];
    }

};

}}

#endif //FILE_CHARGES

