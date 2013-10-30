/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_VSSMGROUP_H_
#define __VOTCA_KMC_VSSMGROUP_H_


#include <votca/kmc/events.h>
#include <votca/kmc/globaleventinfo.h>

namespace votca { namespace kmc {
  
using namespace std;

class Vssmgroup {
    
public:
    
    void Perform_one_step();
    
    votca::tools::Random2 *RandomVariable   
 
    double prob_sum;
    
private:

    
};

/*       double rand_u = 1-RandomVariable->rand_uniform();
        while(rand_u == 0)
        {
            cout << "WARNING: encountered 0 as a random variable! New try." << endl;
            rand_u = 1-RandomVariable->rand_uniform();
        }
        dt = -1 / cumulated_rate * log(rand_u);*/

void Vssmgroup::Compute_sum(Events* events,Globaleventinfo* globevent);


void Vssmgroup::Perform_one_step(Events* events, votca::tools::Random2 *RandomVariable){

    double randn = 1 - RandomVariable->rand_uniform();       
    events->El_non_injection_rates.compute_sum();
    events->Ho_non_injection_rates.compute_sum();
    events->El_injection_rates.compute_sum();
    events->Ho_injection_rates.compute_sum();
}



}} 

#endif
