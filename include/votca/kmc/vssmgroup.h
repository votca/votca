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
    
    void Recompute_in_device(Events* events);
    void Recompute_in_bulk(Events* events);
    double Timestep(votca::tools::Random2 *RandomVariable);
    void Perform_one_step_in_device(Events* events,Graph* graph, State* state, Globaleventinfo* globevent, votca::tools::Random2 *RandomVariable);
    void Perform_one_step_in_bulk(Events* events,Graph* graph, State* state, Globaleventinfo* globevent, votca::tools::Random2 *RandomVariable);
    
    votca::tools::Random2 *RandomVariable;   
 

    
private:

    double tot_probsum;
    
    double non_inject_probsum;
    double inject_probsum;
    
};

double Vssmgroup::Timestep(votca::tools::Random2 *RandomVariable){

    double timestep;
    
    double rand_u = 1-RandomVariable->rand_uniform();
    while(rand_u == 0) {
        rand_u = 1-RandomVariable->rand_uniform();
    }
        
    timestep = -1/tot_probsum*log(rand_u);
    return timestep;
    
}

void Vssmgroup::Recompute_in_device(Events* events){

    non_inject_probsum = events->Non_injection_rates->compute_sum();
    inject_probsum = events->Injection_rates->compute_sum();       
    tot_probsum = non_inject_probsum + inject_probsum;
}

void Vssmgroup::Recompute_in_bulk(Events* events){

    tot_probsum = events->Non_injection_rates->compute_sum();

}

void Vssmgroup::Perform_one_step_in_device(Events* events, Graph* graph, State* state, Globaleventinfo* globevent, votca::tools::Random2 *RandomVariable){

    double randn = RandomVariable->rand_uniform();    
    
    long event_ID;
    Event* chosenevent;
    
    if(randn<inject_probsum/tot_probsum) { // injection event
        randn *= inject_probsum;
        event_ID = events->Injection_rates->search(randn);
        chosenevent = events->Injection_events[event_ID];
    }
    else {
        randn -= inject_probsum/tot_probsum;
        randn *= non_inject_probsum;
        event_ID = events->Non_injection_rates->search(randn);
        chosenevent = events->Non_injection_events[event_ID];        
    }

    events->On_execute(chosenevent, graph, state, globevent);

}

void Vssmgroup::Perform_one_step_in_bulk(Events* events, Graph* graph, State* state, Globaleventinfo* globevent, votca::tools::Random2 *RandomVariable){

    double randn = RandomVariable->rand_uniform();    
    
    long event_ID;
    Event* chosenevent;

    randn *= tot_probsum;
    event_ID = events->Non_injection_rates->search(randn);
    chosenevent = events->Non_injection_events[event_ID];

    events->On_execute(chosenevent, graph, state, globevent);
}


}} 

#endif
