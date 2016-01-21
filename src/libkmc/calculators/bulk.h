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

#ifndef __VOTCA_KMC_BULK_H
#define	__VOTCA_KMC_BULK_H

#include <iostream>

#include <votca/xtp/vssmgroup.h>
#include <votca/xtp/graphkmc.h>
#include <votca/xtp/node.h>
#include <votca/xtp/state.h>
#include <votca/xtp/eventinfo.h>
#include <votca/xtp/event.h>
#include <votca/xtp/events.h>
#include <votca/xtp/profile.h>
#include <votca/xtp/numoutput.h>

using namespace std;

namespace votca { namespace xtp {

    
class Bulk : public KMCCalculator 
{
public:
    
    votca::tools::Random2* randomvariable;
    GraphKMC* graph;
    StateReservoir* state;
    Events* events;
    Vssmgroup* vssmgroup;
    Eventinfo* eventinfo;
    Bsumtree* non_injection_rates;
    Bsumtree* left_injection_rates;
    Bsumtree* right_injection_rates;
    Bsumtree* site_inject_probs;
    Numoutput* numoutput;
    Longrange* longrange;
    
    Bulk() {};
   ~Bulk() {};

    void Initialize(const char *filename, Property *options, const char *outputfile);
    bool EvaluateFrame();

    double sim_time;
    
protected:
   void RunKMC(void); 
            
   
   
};


void Bulk::Initialize(const char *filename, Property *options, const char *outputfile) 
{

    std::cout << "===================================================" << "\n";
    std::cout << "= Initialization phase                            =" << "\n";
    std::cout << "===================================================" << "\n";
    
    eventinfo = new Eventinfo();
    eventinfo->device = false;
    eventinfo->Read_bulk(options);
    
    //Setup random number generator
    srand(eventinfo->random_seed); // srand expects any integer in order to initialise the random number generator
    randomvariable = new votca::tools::Random2();
    randomvariable->init(rand(), rand(), rand(), rand());     
    
    graph = new GraphKMC();
    graph->Initialize_sql(filename);
    std::cout << "number of nodes before graph manipulations: " << graph->Numberofnodes() << "\n";
    std::cout << "simulation box size before graph manipulations: " << graph->Determine_Sim_Box_Size() << "\n";
    
    graph->Setup_bulk_graph(eventinfo);
   
    eventinfo->Graph_Parameters(graph->hopdist(), graph->mindist(), graph->simboxsize(), graph->maxpairdegree(),graph->Average_hole_node_energy(), graph->stddev_hole_node_energy(), graph->Average_electron_node_energy(), graph-> Hole_inject_reorg(), graph->Electron_inject_reorg());
    
    std::cout << "graph object initialized" << "\n";
    std::cout << "max pair degree: " << graph->maxpairdegree() << "\n";
    std::cout << "hopping distance: " << graph->hopdist() << "\n";
    std::cout << "simulation box size: " << graph->simboxsize() << "\n";
    std::cout << "number of nodes " << graph->Numberofnodes() << "\n";    
    
    state = new StateReservoir();
    state->InitStateReservoir();
    std::cout << "state reservoir object initialized" << "\n";
    
    site_inject_probs = new Bsumtree();
    site_inject_probs->initialize(graph->Numberofnodes());
    std::cout << "Random node injector tree initialized" << "\n";
        
    //Random charge distribution
    int nrholes = 0;
    int nrelectrons = 0;
    
    nrholes = ceil(eventinfo->hole_density*graph->Numberofnodes()); 
    nrelectrons = ceil(eventinfo->electron_density*graph->Numberofnodes());
        
    if(nrholes!=0 && nrelectrons != 0) { std::cout << "Double carrier simulation" << "\n";} 
    else { std::cout << "Single carrier simulation" << "\n";} 
    if(nrholes+nrelectrons>graph->Numberofnodes()){ std::cout<< "WARNING: number of electrons and holes is larger than the number of available nodes" << "\n";}
    if(nrholes+nrelectrons == 0) { std::cout << "WARNING: without electrons and/or holes nothing happens";}    
    
    state->Random_init_injection(nrelectrons, nrholes, site_inject_probs, graph, eventinfo, randomvariable);
    std::cout << "randomly placed : " << nrelectrons + nrholes << " charges of which " << nrelectrons << " electrons and " << nrholes << " holes" << "\n";
   
    
    if(state->ReservoirEmpty()) state->Grow(eventinfo->growsize, eventinfo->maxpairdegree);
    std::cout << "state grown after random charge injection" << "\n";
    
    
    non_injection_rates = new Bsumtree();
    left_injection_rates = new Bsumtree();
    right_injection_rates = new Bsumtree();
    std::cout << "rate binary tree initialized" << "\n";

    longrange = new Longrange();    
    
    events = new Events();
    if(eventinfo->coulomb_strength > 0.0) {
        events->Init_non_injection_meshes(eventinfo);
        std::cout << "mesh structure for Coulomb interaction calculations initialized" << "\n";
    }
    events->Initialize_bulk_eventvector(graph,state,eventinfo);
    std::cout << "Initialize (bulk) eventvector (setting number of injection events to zero)" << "\n";
    events->Initialize_rates(non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);
    std::cout << "Fill rate binary trees" << "\n";
    events->Initialize_after_charge_placement(graph,state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);
    std::cout << "Initialize event and rates object after initial placement of charges" << "\n";
    
    if(eventinfo->carrier_trajectory) {state->Init_trajectory(eventinfo->trajectory_filename); std::cout << "Trajectory file initialized" << "\n";}
    
    vssmgroup = new Vssmgroup();
    std::cout << "vssm group initialized" << "\n";

    numoutput = new Numoutput();
    numoutput->Initialize(eventinfo);
    std::cout << "output object initialized" << "\n";


}

bool Bulk::EvaluateFrame() {
    
    // register all QM packages (Gaussian, turbomole, etc))
    // EventFactory::RegisterAll(); 
        
    RunKMC();
   
    delete randomvariable;
    delete state;
    delete events;
    delete graph;
    delete eventinfo;
    delete vssmgroup;
    delete non_injection_rates;
    delete left_injection_rates;
    delete right_injection_rates;
    delete numoutput;
    delete site_inject_probs;
    delete longrange;

    exit(0);
}

void Bulk::RunKMC() {

    
    // to check whether anti repeating methods are useful
    if(eventinfo->repeat_counting) numoutput->Repeat_count_init();
    
    // convergence criteria
    //bool direct_iv_convergence = false;
    //bool direct_reco_convergence = false;
    //int direct_iv_counter = 0; //if the convergence criterium is counted ten times in a row, result is converged
    //int direct_reco_counter = 0;

    if(eventinfo->filament_visualisation)  numoutput->Init_visualisation(graph, eventinfo);
    if(eventinfo->carrier_trajectory) eventinfo->number_of_equilibration_steps = 0;    
    
    std::cout << "total link z distance : " << graph->total_link_distance_z() << "\n";
    std::cout << "average hole site energy : " << eventinfo->avholeenergy << "\n";
    std::cout << "average electron site energy : " << eventinfo->avelectronenergy << "\n"; 
    std::cout << "standard deviation (disorder strength) of hole site energies: " << graph->stddev_hole_node_energy() << "\n";
    std::cout << "standard deviation (disorder strength) of electron site energies: " << graph->stddev_electron_node_energy() << "\n";
    std::cout << "\n";
    std::cout << "===================================================" << "\n";
    std::cout << "= Start of bulk simulation                        =" << "\n";
    std::cout << "===================================================" << "\n";
    std::cout << "\n";
    
    sim_time = 0.0;
    for (long it = 0; it < eventinfo->number_of_equilibration_steps + eventinfo->number_of_steps; it++) {
    
        vssmgroup->Recompute_bulk(non_injection_rates);

        double timestep = vssmgroup->Timestep(randomvariable);
        sim_time += timestep;

        Event* chosenevent;
        chosenevent = vssmgroup->Choose_event_bulk(events, non_injection_rates, randomvariable);
        
        if(eventinfo->filament_visualisation && it <= eventinfo->visualisation_at_nr_steps) numoutput->Update_visualisation(chosenevent);
        if(eventinfo->filament_visualisation && it == eventinfo->visualisation_at_nr_steps) numoutput->Print_visualisation();

        if(eventinfo->carrier_trajectory && ldiv(it,eventinfo->number_of_trajectory_report_steps).rem == 0) state->Print_trajectory(sim_time);

        // check for direct repeats
        if(eventinfo->repeat_counting) numoutput->Repeat_count_update(chosenevent);
      
        if(!eventinfo->carrier_trajectory) numoutput->Update(chosenevent, sim_time, timestep); 
        events->On_execute(chosenevent, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);
        
        // equilibration
   
        if(!eventinfo->carrier_trajectory && it == eventinfo->number_of_equilibration_steps) {
            numoutput->Init_convergence_check(sim_time);
            numoutput->Initialize_equilibrate(eventinfo);
            sim_time = 0.0;
        }
        // convergence checking
        
        if(!eventinfo->carrier_trajectory && (ldiv(it,eventinfo->number_of_report_steps ).rem==0 && it> eventinfo->number_of_equilibration_steps)) numoutput->Convergence_check(sim_time, eventinfo);

        // direct output
        if(!eventinfo->carrier_trajectory && ldiv(it,eventinfo->number_of_report_steps ).rem==0){
            numoutput->Write(it, sim_time, timestep, eventinfo);
            std::cout << "\n";
        }
        
        if(eventinfo->carrier_trajectory && ldiv(it,eventinfo->number_of_report_steps ).rem==0){
            std::cout << "step " << it << " of " << eventinfo->number_of_equilibration_steps + eventinfo->number_of_steps << " timesteps" << "\n";
        }        
        
        // break out of loop
        if(!eventinfo->carrier_trajectory && numoutput->iv_conv()) {break;}
        
    }

}

}}


#endif	/* __VOTCA_KMC_BULK_H */
