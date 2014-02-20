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
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/event.h>
#include <votca/kmc/events.h>
#include <votca/kmc/profile.h>
#include <votca/kmc/numoutput.h>

using namespace std;

namespace votca { namespace kmc {
    
//typedef votca::tools::vec myvec;

// LUMO level substracted!!!
   
class Diode : public KMCCalculator 
{
public:
    
    votca::tools::Random2 *RandomVariable;
    GraphDevice* graph;
    StateDevice* state;
    Events* events;
    Vssmgroup* vssmgroup;
    Eventinfo* eventdata;
    Longrange* longrange;
    Bsumtree* non_injection_rates;
    Bsumtree* injection_rates;
    Bsumtree* site_inject_probs;
    Numoutput* numoutput;
    
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

    //Setup random number generator
    srand(eventdata->seed); // srand expects any integer in order to initialise the random number generator
    RandomVariable = new votca::tools::Random2();
    RandomVariable->init(rand(), rand(), rand(), rand());     
    
    cout << "Initializing" << endl;
    
    eventdata = new Eventinfo();
    eventdata->Read(options);
    
    graph = new GraphDevice();
    graph->Initialize(filename);
    graph->Setup_device_graph(eventdata->left_electrode_distance, eventdata->right_electrode_distance, true, eventdata);
    eventdata->Graph_Parameters(graph->hopdist(), graph->simboxsize(), graph->maxpairdegree(),graph->Av_hole_node_energy());
    eventdata->Set_field(); // convert voltage to electric field

    std::cout << "graph initialized" << endl;
    std::cout << "max pair degree: " << graph->maxpairdegree() << endl;
    std::cout << "hopping distance: " << graph->hopdist() << endl;
    std::cout << "simulation box size: " << graph->simboxsize() << endl;
    std::cout << "number of left electrode injector nodes " << graph->left()->links().size() << endl;
    std::cout << "number of right electrode injector nodes " << graph->right()->links().size() << endl;    
    
    longrange = new Longrange(graph,eventdata);
    longrange->Initialize(eventdata);
    
    std::cout << "longrange profile initialized" << endl;

    state = new StateDevice();
    state->InitStateDevice();
    
    site_inject_probs = new Bsumtree();
    site_inject_probs->initialize(graph->Numberofnodes()); // take care of electrode nodes
    state->Random_init_injection((int) Hole, site_inject_probs, graph, eventdata, RandomVariable);
    
    if(state->ReservoirEmpty()) state->Grow(eventdata->growsize, eventdata->maxpairdegree);
    
    std::cout << "state initialized" << endl;
    
    non_injection_rates = new Bsumtree();
    injection_rates = new Bsumtree();
    
    std::cout << "binary tree structures initialized" << endl;

    events = new Events();
    events->Init_non_injection_meshes(eventdata);
    events->Initialize_eventvector(graph,state,longrange,eventdata);
    events->Initialize_rates(non_injection_rates,injection_rates,eventdata);
    events->Init_injection_meshes(state, eventdata);
    events->Initialize_after_charge_placement(graph,state, longrange, non_injection_rates, injection_rates, eventdata);
    std::cout << "event vectors and meshes initialized" << endl;

    vssmgroup = new Vssmgroup();
    
    numoutput = new Numoutput();
    numoutput->Initialize();

}

bool Diode::EvaluateFrame() {
    
    // register all QM packages (Gaussian, turbomole, etc))
    // EventFactory::RegisterAll(); 
        
    RunKMC();
    
    delete RandomVariable;
    delete state;
    delete events;
    delete graph;
    delete eventdata;
    delete longrange;
    delete vssmgroup;
    delete non_injection_rates;
    delete injection_rates;
    delete numoutput;
    delete site_inject_probs;
    exit(0);
}

void Diode::RunKMC() {
    
    // to check whether anti repeating methods are useful
    int repeat_counter = 0; 
    int old_from_node_id = -10;
    int old_to_node_id = -10;
    
    // convergence criteria
    bool direct_iv_convergence = false;
    bool direct_reco_convergence = false;
    int direct_iv_counter = 0; //if the convergence criterium is counted ten times in a row, result is converged
    int direct_reco_counter = 0;
    
    sim_time = 0.0;
    for (long it = 0; it < 2*eventdata->nr_equilsteps + eventdata->nr_timesteps; it++) {
        // Update longrange cache (expensive, so not done at every timestep)
        if(ldiv(it, eventdata->steps_update_longrange).rem == 0 && it>0){
            longrange->Update_cache(eventdata);
            events->Recompute_all_events(state, longrange, non_injection_rates, injection_rates, eventdata);
        }
       
        vssmgroup->Recompute_device(non_injection_rates, injection_rates);
        double timestep = vssmgroup->Timestep(RandomVariable);
        sim_time += timestep;
 
        Event* chosenevent = vssmgroup->Choose_event_device(events, non_injection_rates, injection_rates, RandomVariable);
        numoutput->Update(chosenevent, sim_time, timestep); 

        events->On_execute(chosenevent, graph, state, longrange, non_injection_rates, injection_rates, eventdata, it);

        // check for direct repeats
        
        int goto_node_id = chosenevent->link()->node2()->id();
        int from_node_id = chosenevent->link()->node1()->id();
        if(goto_node_id == old_from_node_id && from_node_id == old_to_node_id) repeat_counter++;
        old_from_node_id = from_node_id;
        old_to_node_id = goto_node_id;

        if(it == 2*eventdata->nr_equilsteps) numoutput->Init_convergence_check(sim_time);
        
        // equilibration
        
        if(it == eventdata->nr_equilsteps || it == 2*eventdata->nr_equilsteps) {
            numoutput->Initialize_equilibrate();
            sim_time = 0.0;
        }
        
        // convergence checking
        
        if(ldiv(it,100).rem==0 && it> 2*eventdata->nr_equilsteps) numoutput->Convergence_check(sim_time, eventdata);

        // direct output
        if(ldiv(it,100).rem==0){
            std::cout << it << " " << repeat_counter << " " << 
                         numoutput->iv_conv() << " " << numoutput->iv_count() << " " << 
                         numoutput->reco_conv() << " " << numoutput->reco_count() <<  " " << 
                         sim_time << " " << timestep << " ";
            numoutput->Write(sim_time);
        }
        
        // break out of loop
        if(numoutput->iv_conv() && numoutput->reco_conv()) {break;}
        
    }

}

}}


#endif	/* __VOTCA_KMC_DIODE_H */