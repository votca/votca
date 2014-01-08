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

#ifndef __VOTCA_KMC_DIODE_H
#define	__VOTCA_KMC_DIODE_H

#include <iostream>

#include <votca/kmc/vssmgroup.h>
#include <votca/kmc/graphdevice.h>
#include <votca/kmc/node.h>
#include <votca/kmc/state.h>
#include <votca/kmc/mesh.h>
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/event.h>
#include <votca/kmc/events.h>
#include <votca/kmc/profile.h>

using namespace std;

namespace votca { namespace kmc {
    
//typedef votca::tools::vec myvec;

   
class Diode : public KMCCalculator 
{
public:
    
    GraphDevice* graph;
    StateDevice* state;
    Events* events;
    Vssmgroup* vssmgroup;
    Eventinfo* eventdata;
    Longrange* longrange;
    Bsumtree* non_injection_rates;
    Bsumtree* injection_rates;
    
    Diode() {};
   ~Diode() {};

    void Initialize(const char *filename, Property *options, const char *outputfile);
    bool EvaluateFrame();

    double sim_time;
    
protected:
   void RunKMC(void); 
            
private:
    static const double kB   = 8.617332478E-5; // eV/K
    static const double hbar = 6.5821192815E-16; // eV*s
    static const double eps0 = 8.85418781762E-12/1.602176565E-19; // e**2/eV/m = 8.85418781762E-12 As/Vm
    static const double epsr = 3.0; // relative material permittivity
    static const double Pi   = 3.14159265358979323846;
   
};


void Diode::Initialize(const char *filename, Property *options, const char *outputfile) {
    
    cout << "Initializing" << endl;
    
    eventdata = new Eventinfo();
    eventdata->Read(options);
    
    graph = new GraphDevice();
    graph->Initialize(filename);
    graph->Setup_device_graph(eventdata->left_electrode_distance, eventdata->right_electrode_distance, eventdata);
    std::cout << "max pair degree: " << graph->maxpairdegree() << endl;
    std::cout << "hopping distance: " << graph->hopdist() << endl;
    std::cout << "simulation box size: " << graph->simboxsize() << endl;
    std::cout << "number of left electrode injector nodes " << graph->left()->links().size() << endl;
    std::cout << "number of right electrode injector nodes " << graph->right()->links().size() << endl;
    eventdata->Graph_Parameters(graph->hopdist(), graph->simboxsize(), graph->maxpairdegree());

    longrange = new Longrange(graph,eventdata);
    longrange->Initialize(eventdata);

    state = new StateDevice();
    state->InitStateDevice();
    state->Grow(eventdata->growsize, eventdata->maxpairdegree); //initial growth

/*    state->Grow(10, graph->maxpairdegree());
    int carrier_ID = state->Buy();
    Carrier* newcarrier = state->GetCarrier(carrier_ID);
    Node* carrier_node = graph->GetNode(20);
    newcarrier->SetCarrierNode(carrier_node);
    newcarrier->SetCarrierType(2);
    carrier_node->AddCarrier(carrier_ID);
    longrange->Init_Load_State(state, eventdata);*/
    
    non_injection_rates = new Bsumtree();
    injection_rates = new Bsumtree();
    events = new Events();
    events->Initialize_eventvector(graph,state,longrange,non_injection_rates,injection_rates,eventdata);

    vssmgroup = new Vssmgroup();

}

bool Diode::EvaluateFrame() {
    
    // register all QM packages (Gaussian, turbomole, etc))
    // EventFactory::RegisterAll(); 
        
    RunKMC();
    
    delete state;
    delete events;
    delete graph;
    delete eventdata;
    delete longrange;
    delete vssmgroup;
    delete non_injection_rates;
    delete injection_rates;    
    exit(0);
}

void Diode::RunKMC() {
    
    //Setup random number generator
    srand(eventdata->seed); // srand expects any integer in order to initialise the random number generator
    votca::tools::Random2 *RandomVariable = new votca::tools::Random2();
    RandomVariable->init(rand(), rand(), rand(), rand());    

    sim_time = 0.0;
    for (long it = 0; it < 2*eventdata->nr_equilsteps + eventdata->nr_timesteps; it++) {
        
        // Update longrange cache (expensive, so not done at every timestep)
        if(ldiv(it, eventdata->steps_update_longrange).rem == 0 && it>0){
            longrange->Update_cache(eventdata);
            events->Recompute_all_events(state, longrange, non_injection_rates, injection_rates, eventdata);
        }
        
        vssmgroup->Recompute_in_device(events, non_injection_rates, injection_rates);
        sim_time += vssmgroup->Timestep(RandomVariable);
        Event* chosenevent = vssmgroup->Choose_event(events, non_injection_rates, injection_rates, RandomVariable);

        events->On_execute(chosenevent, graph, state, longrange, non_injection_rates, injection_rates, eventdata);

        // take care of output here

    }
}

}}


#endif	/* __VOTCA_KMC_DIODE_H */