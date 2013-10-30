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

    double el_probsum;
    double ho_probsum;
    double tot_probsum;
    
    double el_non_inject_probsum;
    double ho_non_inject_probsum;
    double el_inject_probsum;
    double ho_inject_probsum;
    
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
    
    if(events->el_dirty) {
        el_non_inject_probsum = events->El_non_injection_rates->compute_sum();
        el_inject_probsum = events->El_injection_rates->compute_sum();       
        el_probsum = el_non_inject_probsum + el_inject_probsum;
        events->el_dirty = false;
    }

    if(events->ho_dirty) {
        ho_non_inject_probsum = events->Ho_non_injection_rates->compute_sum();
        ho_inject_probsum = events->Ho_injection_rates->compute_sum();       
        ho_probsum = ho_non_inject_probsum + ho_inject_probsum;
        events->ho_dirty = false;
    }

    tot_probsum = el_probsum + ho_probsum;
}

void Vssmgroup::Recompute_in_bulk(Events* events){
    
    if(events->el_dirty) {
        el_probsum = events->El_non_injection_rates->compute_sum();
        events->el_dirty = false;
    }

    if(events->ho_dirty) {
        ho_probsum = events->Ho_non_injection_rates->compute_sum();
        events->ho_dirty = false;
    }

    tot_probsum = el_probsum + ho_probsum;
    
}

void Vssmgroup::Perform_one_step_in_device(Events* events, Graph* graph, State* state, Globaleventinfo* globevent, votca::tools::Random2 *RandomVariable){

    double randn = RandomVariable->rand_uniform();    
    
    long event_ID;
    Event* chosenevent;
    
    if(randn<el_probsum/tot_probsum) { // electron event
        if(randn< el_non_inject_probsum/tot_probsum) { // non injection event
            randn *= el_non_inject_probsum;
            event_ID = events->El_non_injection_rates->search(randn);
            chosenevent = events->El_non_injection_events[event_ID];
        }
        else { // injection event
            randn -= el_non_inject_probsum/tot_probsum;
            randn *= el_inject_probsum;
            event_ID = events->El_injection_rates->search(randn);
            chosenevent = events->El_injection_events[event_ID];
        }
    }
    else { //Hole event
        randn -= el_probsum/tot_probsum;
        if(randn< ho_non_inject_probsum/tot_probsum) { // non injection event
            randn *= ho_non_inject_probsum;
            event_ID = events->Ho_non_injection_rates->search(randn);
            chosenevent = events->Ho_non_injection_events[event_ID];
        }
        else { // injection event
            randn -= ho_non_inject_probsum/tot_probsum;
            randn *= ho_inject_probsum;
            event_ID = events->Ho_injection_rates->search(randn);
            chosenevent = events->Ho_injection_events[event_ID];
        }
    }
    events->On_execute(chosenevent, graph, state, globevent);
}

void Vssmgroup::Perform_one_step_in_bulk(Events* events, Graph* graph, State* state, Globaleventinfo* globevent, votca::tools::Random2 *RandomVariable){

    double randn = RandomVariable->rand_uniform();    
    
    //first search the particular bsumtree
   
    long event_ID;
    Event* chosenevent;
    
    if(randn<el_probsum/tot_probsum) { // electron event
        randn *= el_probsum;
        event_ID = events->El_non_injection_rates->search(randn);
        chosenevent = events->El_non_injection_events[event_ID];
    }
    else { //Hole event
        randn -= el_probsum/tot_probsum;
        randn *= ho_probsum;
        event_ID = events->Ho_non_injection_rates->search(randn);
        chosenevent = events->Ho_non_injection_events[event_ID];
    }
    events->On_execute(chosenevent, graph, state, globevent);
}


}} 

#endif
