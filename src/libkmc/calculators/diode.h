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
//#include <votca/kmc/graphsql.h>
//#include <votca/kmc/graphcubic.h>

//#include <votca/kmc/vssmgroup.h>
#include <votca/kmc/graphdevice.h>
#include <votca/kmc/node.h>
#include <votca/kmc/state.h>
#include <votca/kmc/mesh.h>
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/event.h>
#include <votca/kmc/events.h>

using namespace std;

namespace votca { namespace kmc {
    
//typedef votca::tools::vec myvec;

   
class Diode : public KMCCalculator 
{
public:
    
    GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph;
    StateDevice* state;
    Events* events;
//    Vssmgroup* vssmgroup;
    Eventinfo* eventdata;
    
    Diode() {};
   ~Diode() {};

    void Initialize(const char *filename, Property *options, const char *outputfile);
    bool EvaluateFrame();

    double sim_time;
    
    // input parameters (put in globaleventinfo?)
//    int seed; long nr_equilsteps; long nr_timesteps; long steps_update_longrange;
//   int nx; int ny; int nz; double lattice_constant; double hopdist; double disorder_strength; 
//    double disorder_ratio; CorrelationType correlation_type; double left_electrode_distance; double right_electrode_distance;
    
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
    
    graph = new GraphDevice<GraphSQL, NodeSQL, LinkSQL>();
    graph->Initialize(filename);
    graph->Setup_device_graph(eventdata->left_electrode_distance, eventdata->right_electrode_distance);
    std::cout << "max pair degree: " << graph->maxpairdegree() << endl;
    std::cout << "hopping distance: " << graph->hopdist() << endl;
    std::cout << "simulation box size: " << graph->simboxsize() << endl;
    std::cout << "number of left electrode injector nodes " << graph->left()->links().size() << endl;
    std::cout << "number of right electrode injector nodes " << graph->right()->links().size() << endl;

    state = new StateDevice();
    state->InitState();
    
    events = new Events();
    
    std::cout << graph->GetNode(10)->occ() << endl;
    state->Grow(10);
    int carrier_ID = state->Buy();
    Carrier* newcarrier = state->GetCarrier(carrier_ID);
    Node* carrier_node = graph->GetNode(20);
    newcarrier->SetCarrierNode(carrier_node);
    newcarrier->SetCarrierType(2);
    carrier_node->AddCarrier(carrier_ID);
    vector<Link*> links = carrier_node->links();
    
    
    Carrier* newcarrier2 = state->GetCarrier(3);
    newcarrier2->SetCarrierType(2);
    
    events->Initialize_eventvector(graph,state,eventdata);

    
    
    
    
/*    for(int it = 0; it < graph->maxpairdegree(); it++) {
        Event* newevent = new Event();
//        std::cout << it << " " << newevent->rate() << " " << newevent->init_type() << " " << newevent->final_type() << " " << endl;
    }
    typename std::vector<Link*>::iterator it;
    for(it = links.begin(); it != links.end(); it++) {
        Event* newevent = new Event((*it), newcarrier->type(), eventdata,state);
//        std::cout << (*it)->id() << " " << newevent->rate() << " " << newevent->init_type() << " " << newevent->final_type() << " " 
//                  << abs((*it)->r12()) << " " << (*it)->r12() << " " <<  exp(-1.0*eventdata->alpha*abs((*it)->r12())) << " " 
//                  << endl;

    }
    vector<Link*> inj_links = graph->left()->links();
    typename std::vector<Link*>:: iterator inj;
    for(inj = inj_links.begin(); inj != inj_links.end(); inj++) {
        Event* newevent = new Event((*inj), 2, eventdata,state);
//        std::cout << (*inj)->id() << " " << newevent->rate() << " " << newevent->init_type() << " " << newevent->final_type() << " " 
//                  << abs((*inj)->r12()) << " " << (*inj)->r12() << " " <<  exp(-1.0*eventdata->alpha*abs((*inj)->r12())) << " " 
//                  << endl;

    }*/
    
    delete state;
    delete events;
    delete graph;
    delete eventdata;    
    exit(0);
    

//    events = new Events();
//    vssmgroup = new Vssmgroup();
}

bool Diode::EvaluateFrame() {
    
    // register all QM packages (Gaussian, turbomole, etc))
    // EventFactory::RegisterAll(); 
        
    RunKMC();
}

void Diode::RunKMC() {
 
    exit(0);
    /*    
    //Setup random number generator
    seed = 1;
    srand(seed); // srand expects any integer in order to initialise the random number generator
    votca::tools::Random2 *RandomVariable = new votca::tools::Random2();
    RandomVariable->init(rand(), rand(), rand(), rand());    
    
    //Initialize all structures
    graph->hopdist = hopdist;
    graph->Generate_cubic_graph(nx, ny, nz, lattice_constant, disorder_strength,RandomVariable, disorder_ratio, 
                                correlation_type, left_electrode_distance, right_electro_distance,globevent);   
    state->Init();    
    events->Initialize_eventvector(graph, state, globevent);
    events->Initialize_longrange (graph, globevent);
    events->Recompute_all_injection_events(graph, globevent);
    events->Recompute_all_non_injection_events(graph, state, globevent); 
    
    
    sim_time = 0.0;
    for (long it = 0; it < 2*nr_equilsteps + nr_timesteps; it++) {
        
        // Update longrange cache (expensive, so not done at every timestep)
        if(ldiv(it, steps_update_longrange).rem == 0 && it>0){
            events->longrange->Update_cache(graph->sim_box_size, globevent);
            events->Recompute_all_injection_events(graph, globevent);
            events->Recompute_all_non_injection_events(graph,state,globevent);
        }
        
        vssmgroup->Recompute_in_device(events);
        sim_time += vssmgroup->Timestep(RandomVariable);
        vssmgroup->Perform_one_step_in_device(events,graph,state,globevent,RandomVariable);
    }
     */
}

}}


#endif	/* __VOTCA_KMC_DIODE_H */