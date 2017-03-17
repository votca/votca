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

#include "votca/xtp/vssmgroup.h"


using namespace std;

namespace votca {
    namespace xtp {
double Vssmgroup::Timestep(votca::tools::Random2 *RandomVariable){
    double timestep;
    
    double rand_u = 1-RandomVariable->rand_uniform();
    while(rand_u == 0) {
        rand_u = 1-RandomVariable->rand_uniform();
    }
       
    timestep = (-1.0/tot_probsum)*log(rand_u);
    return timestep;
    
}

void Vssmgroup::Recompute_device(Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates){

    non_inject_probsum = non_injection_rates->compute_sum();
    if(left_injection_rates->getnrrates() != 0) {
        left_inject_probsum = left_injection_rates->compute_sum();
    }
    else {
        left_inject_probsum = 0.0;
    }
    if(right_injection_rates->getnrrates() != 0) {
        right_inject_probsum = right_injection_rates->compute_sum();
    }
    else {
        right_inject_probsum = 0.0;
    }
    tot_probsum = non_inject_probsum + left_inject_probsum + right_inject_probsum;
//   std::cout << "non " << non_inject_probsum << " left " << left_inject_probsum << " right " << right_inject_probsum << endl;
}

void Vssmgroup::Recompute_bulk(Bsumtree* non_injection_rates){

    tot_probsum = non_injection_rates->compute_sum();
}

void Vssmgroup::Recompute_injection(Bsumtree* injection_rates){
    
    tot_probsum = injection_rates->compute_sum();
}

Event* Vssmgroup::Choose_event_device(Events* events, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, votca::tools::Random2 *RandomVariable){

    double randn = 0.0;
    while(randn == 0.0) {
        randn = tot_probsum*RandomVariable->rand_uniform();
    }
    
    long event_ID;
    Event* chosenevent;
    if(randn<left_inject_probsum) { // injection event
        event_ID = left_injection_rates->search(randn);
        chosenevent = events->get_injection_event(0,event_ID);
    }
    else {
        randn -= left_inject_probsum;
        if(randn<right_inject_probsum){
            event_ID = right_injection_rates->search(randn);
            chosenevent = events->get_injection_event(1,event_ID);       
        }
        else {   
            randn -= right_inject_probsum;
            event_ID = non_injection_rates->search(randn);
            chosenevent = events->get_non_injection_event(event_ID);       
        }
    }
    return chosenevent;
}

Event* Vssmgroup::Choose_event_bulk(Events* events, Bsumtree* non_injection_rates,votca::tools::Random2 *RandomVariable){

    double randn = 0.0;
    while(randn == 0.0) {
        randn = tot_probsum*RandomVariable->rand_uniform();
    }
    
    long event_ID;

    Event* chosenevent;
    event_ID = non_injection_rates->search(randn);
    chosenevent = events->get_non_injection_event(event_ID);       
    return chosenevent;

}

Event* Vssmgroup::Choose_injection_event(Events* events,  int electrodeID, Bsumtree* injection_rates,votca::tools::Random2 *RandomVariable){

    double randn = 0.0;
    while(randn == 0.0) {
        randn = tot_probsum*RandomVariable->rand_uniform();
    }
    long event_ID;
    Event* chosenevent;
    event_ID = injection_rates->search(randn);
    chosenevent = events->get_injection_event(electrodeID, event_ID);       
    return chosenevent;

}


        
    }
}
