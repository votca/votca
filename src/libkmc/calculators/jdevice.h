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

#ifndef __VOTCA_KMC_JDEVICE_H
#define	__VOTCA_KMC_JDEVICE_H

#include <iostream>

#include <votca/kmc/vssmgroup.h>
#include <votca/kmc/graphkmc.h>
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
   
class Jdevice : public KMCCalculator 
{
public:
    
    votca::tools::Random2 *RandomVariable;
    GraphKMC* graph;
    StateReservoir* state;
    Events* events;
    Vssmgroup* vssmgroup;
    Eventinfo* eventinfo;
    Longrange* longrange;
    Bsumtree* non_injection_rates;
    Bsumtree* left_injection_rates;
    Bsumtree* right_injection_rates;
    Bsumtree* site_inject_probs;

    Numoutput* numoutput;
    
    Jdevice() {};
   ~Jdevice() {};

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


void Jdevice::Initialize(const char *filename, Property *options, const char *outputfile) {

    //Setup random number generator
    srand(eventinfo->seed); // srand expects any integer in order to initialise the random number generator
    RandomVariable = new votca::tools::Random2();
    RandomVariable->init(rand(), rand(), rand(), rand());     
    
    cout << "Initializing" << endl;
    
    eventinfo = new Eventinfo();
    eventinfo->Read_device(options);
    
    graph = new GraphKMC();
    graph->Initialize(filename);
    std::cout << "number of nodes: " << graph->Numberofnodes() << endl;
    if(graph->el_reorg()) { std::cout << "WARNING: zero electron reorganization energy" << endl;}
    if(graph->ho_reorg()) { std::cout << "WARNING: zero hole reorganization energy" << endl;}
    
    
    if(eventinfo->device > 0){
        graph->Setup_device_graph(eventinfo->left_electrode_distance, eventinfo->right_electrode_distance, eventinfo->resize, eventinfo);
    }
    else {
        graph->Setup_bulk_graph(eventinfo->resize, eventinfo);
    }    
   
    eventinfo->Graph_Parameters(graph->hopdist(), graph->mindist(), graph->simboxsize(), graph->maxpairdegree(),graph->Average_hole_node_energy(), graph->Average_electron_node_energy());
//    eventinfo->Set_field(eventinfo->voltage* (eventinfo->simboxsize.x()/(eventinfo->simboxsize.x() + eventinfo->left_oxide_thickness + eventinfo->right_oxide_thickness))); // convert voltage to electric field
    if(eventinfo->device > 0) eventinfo->Set_field(eventinfo->voltage); // convert voltage to electric field

    
    std::cout << "graph initialized" << endl;
    std::cout << "max pair degree: " << graph->maxpairdegree() << endl;
    std::cout << "hopping distance: " << graph->hopdist() << endl;
    std::cout << "simulation box size: " << graph->simboxsize() << endl;
    if(eventinfo->device > 0) {
        std::cout << "number of left electrode injector nodes " << graph->left()->links().size() << endl;
        std::cout << "number of right electrode injector nodes " << graph->right()->links().size() << endl;
    }
    std::cout << "number of nodes " << graph->Numberofnodes() << endl;    
    
    if(eventinfo->device >0) {    
        longrange = new Longrange(graph,eventinfo);
        if(eventinfo->longrange_slab) longrange->Initialize_slab(graph,eventinfo);
        else longrange->Initialize(eventinfo);
    }
    else {
        longrange = new Longrange();
    }
    
    std::cout << "longrange profile initialized" << endl;

    state = new StateReservoir();
    state->InitStateReservoir();
    
    site_inject_probs = new Bsumtree();
    site_inject_probs->initialize(graph->Numberofnodes()); // take care of electrode nodes
    
    //Random charge distribution (assume homogeneous distributed over device)
    int nrcharges;
    if(eventinfo->device>0) nrcharges = (eventinfo->voltage)/(2*Pi*eventinfo->coulomb_strength*eventinfo->simboxsize.x()*eventinfo->simboxsize.x())*graph->Numberofnodes();
    if(eventinfo->device==0) nrcharges = eventinfo->ho_density*graph->Numberofnodes();
    if(eventinfo->traj_store) {
        int nrelectrons = eventinfo->nr_electrons;
        int nrholes = eventinfo->nr_holes;
        if(nrholes != 0) {
            nrcharges = nrholes;
            state->Random_init_injection((int) Hole, nrcharges, site_inject_probs, graph, eventinfo, RandomVariable);
        }
        else if(nrelectrons != 0) {
            nrcharges = nrelectrons;
            state->Random_init_injection((int) Electron, nrcharges, site_inject_probs, graph, eventinfo, RandomVariable);
        }
    }
    else {
        state->Random_init_injection((int) Hole, nrcharges, site_inject_probs, graph, eventinfo, RandomVariable);
    }
    std::cout << "initial nrcharges: " << nrcharges << endl;
    
    std::cout << "charges injected" << endl;
    
    if(state->ReservoirEmpty()) state->Grow(eventinfo->growsize, eventinfo->maxpairdegree);
    
    std::cout << "state initialized" << endl;
    
    non_injection_rates = new Bsumtree();
    left_injection_rates = new Bsumtree();
    right_injection_rates = new Bsumtree();

    std::cout << "binary tree structures initialized" << endl;
    
    events = new Events();
    events->Init_non_injection_meshes(eventinfo);
    events->Initialize_device_eventvector(graph,state,longrange,eventinfo);
    events->Initialize_rates(non_injection_rates, left_injection_rates, right_injection_rates,eventinfo);
    if(eventinfo->device>0) events->Init_injection_meshes(state, eventinfo);
    events->Initialize_after_charge_placement(graph,state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);
    
    std::cout << "event vectors and meshes initialized" << endl;

    vssmgroup = new Vssmgroup();

    numoutput = new Numoutput();
    numoutput->Initialize();

}

bool Jdevice::EvaluateFrame() {
    
    // register all QM packages (Gaussian, turbomole, etc))
    // EventFactory::RegisterAll(); 
        
    RunKMC();
    delete RandomVariable;
    delete state;
    delete events;
    delete graph;
    delete eventinfo;
    delete longrange;
    delete vssmgroup;
    delete non_injection_rates;
    delete left_injection_rates;
    delete right_injection_rates;
    delete numoutput;
    delete site_inject_probs;
    exit(0);
}

void Jdevice::RunKMC() {
    
    // to check whether anti repeating methods are useful
    int repeat_counter = 0; 
    int old_from_node_id = -10;
    int old_to_node_id = -10;
    
    // convergence criteria
    bool direct_iv_convergence = false;
    bool direct_reco_convergence = false;
    int direct_iv_counter = 0; //if the convergence criterium is counted ten times in a row, result is converged
    int direct_reco_counter = 0;

    int numcharges_distrib;
    
    if(eventinfo->device >0) {
        for(int it= 0; it< eventinfo->number_layers; it++ ){
            std::cout << longrange->number_of_nodes(it) << " ";
        }
        std::cout << endl;
    }
    
    std::cout << "total link x distance : " << graph->total_link_distance_x() << endl;
    std::cout << "average hole energy : " << eventinfo->avholeenergy << endl;
    std::cout << "disorder strength: " << graph->stddev_hole_node_energy() << endl;
    
    if(eventinfo->viz_store)  numoutput->Init_visualisation(graph, eventinfo);

    
    sim_time = 0.0;
    for (long it = 0; it < 2*eventinfo->nr_equilsteps + eventinfo->nr_timesteps; it++) {
        //    for (long it = 0; it < 100; it++) {
        
        if(eventinfo->device>0) {
            if(ldiv(it, eventinfo->steps_update_longrange).rem == 0 && it>0){
                if(eventinfo->longrange_slab) longrange->Update_cache_slab(graph,eventinfo);
                else                          longrange->Update_cache(eventinfo);

    //            double left_oxide_drop = longrange->Calculate_left_oxide_layer(eventinfo);
    //            double right_oxide_drop = longrange->Calculate_right_oxide_layer(eventinfo);
    //            double in_organic_voltage = eventinfo->voltage - left_oxide_drop - right_oxide_drop;
    //            std::cout << "APPLIED VOLTAGE: " << eventinfo->voltage << ", VOLTAGE AFTER OXIDE DROP " << in_organic_voltage << endl;

    //            eventinfo->Set_field(in_organic_voltage);

                // update extrapolation of potential drop


                events->Recompute_all_events(state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);
            }
        }
        if(eventinfo->device == 0) vssmgroup->Recompute_bulk(non_injection_rates);
        if(eventinfo->device == 1) vssmgroup->Recompute_device(non_injection_rates, left_injection_rates, right_injection_rates);

        double timestep = vssmgroup->Timestep(RandomVariable);
        sim_time += timestep;
        
        Event* chosenevent;
        if(eventinfo->device == 0) chosenevent = vssmgroup->Choose_event_bulk(events, non_injection_rates, RandomVariable);
        if(eventinfo->device == 1) chosenevent = vssmgroup->Choose_event_device(events, non_injection_rates, left_injection_rates, right_injection_rates, RandomVariable);
        
        if(eventinfo->viz_store && it <= eventinfo->viz_nr_timesteps) numoutput->Update_visualisation(chosenevent);
        if(eventinfo->viz_store && it == eventinfo->viz_nr_timesteps) numoutput->Print_visualisation();
        
        numoutput->Update(chosenevent, sim_time, timestep); 
        events->On_execute(chosenevent, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);

        
        // check for direct repeats        
        int goto_node_id = chosenevent->link()->node2()->id();
        int from_node_id = chosenevent->link()->node1()->id();
        if(goto_node_id == old_from_node_id && from_node_id == old_to_node_id) repeat_counter++;
        old_from_node_id = from_node_id;
        old_to_node_id = goto_node_id;

        
         if(it == eventinfo->nr_equilsteps || it == 2*eventinfo->nr_equilsteps) numoutput->Init_convergence_check(sim_time);
        
        // equilibration
   
        if(it == eventinfo->nr_equilsteps || it == 2*eventinfo->nr_equilsteps) {
            numoutput->Initialize_equilibrate();
            sim_time = 0.0;
        }
        // convergence checking
        
        if(ldiv(it,10000).rem==0 && it> 2*eventinfo->nr_equilsteps) numoutput->Convergence_check(sim_time, eventinfo);

        // direct output
        if(ldiv(it,10000).rem==0){
            std::cout << it << " " << repeat_counter << " " << 
                         numoutput->iv_conv() << " " << numoutput->iv_count() << " " << 
                         numoutput->reco_conv() << " " << numoutput->reco_count() <<  " " << 
                         sim_time << " " << timestep << " " << vssmgroup->totprobsum() << " "  << vssmgroup->noninjectprobsum() << " "  << vssmgroup->leftinjectprobsum() << " "  << vssmgroup->rightinjectprobsum() << " " ;
            numoutput->Write(sim_time);
            std::cout << endl;
//            std::cout << "position" << endl;
//           
//            for(int i =0; i< eventinfo->number_layers; i++) {
//                std::cout << longrange->position(i) << " ";
//            }
//            
//            std::cout << endl;
            if(eventinfo->device>0) {
                std::cout << "charge" << endl;
           
                for(int i =0; i< eventinfo->number_layers; i++) {
                    std::cout << longrange->Get_cached_density(i, eventinfo) << " ";
                }
            
                std::cout << endl;
            }
//            std::cout << "pot" << endl;
//           
//            for(int i =0; i< eventinfo->number_layers; i++) {
//                std::cout << longrange->Get_cached_longrange(i) << " ";
//            }            
//            std::cout << endl;
        }
        
        // break out of loop
//        if(numoutput->iv_conv() && numoutput->reco_conv()) {break;}
        if(numoutput->iv_conv()) {break;}
        
    }

}

}}


#endif	/* __VOTCA_KMC_JDEVICE_H */
