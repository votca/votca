/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef FILE_CHARGES
#define FILE_CHARGES

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>

namespace votca { namespace xtp {

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

