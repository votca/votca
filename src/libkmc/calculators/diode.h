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
    Vssmgroup* left_chargegroup;
    Vssmgroup* right_chargegroup;
    Eventinfo* eventdata;
    Longrange* longrange;
    Bsumtree* non_injection_rates;
    Bsumtree* left_injection_rates;
    Bsumtree* right_injection_rates;
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
    std::cout << "number of nodes: " << graph->Numberofnodes() << endl;
    if(graph->el_reorg()) { std::cout << "WARNING: zero electron reorganization energy" << endl;}
    if(graph->ho_reorg()) { std::cout << "WARNING: zero hole reorganization energy" << endl;}
    
    
    if(eventdata->device > 0){
        graph->Setup_device_graph(eventdata->left_electrode_distance, eventdata->right_electrode_distance, eventdata->resize, eventdata);
    }
    else {
        graph->Setup_bulk_graph(eventdata->resize, eventdata);
    }    
   
    eventdata->Graph_Parameters(graph->hopdist(), graph->mindist(), graph->simboxsize(), graph->maxpairdegree(),graph->Average_hole_node_energy(), graph->Average_electron_node_energy());
//    eventdata->Set_field(eventdata->voltage* (eventdata->simboxsize.x()/(eventdata->simboxsize.x() + eventdata->left_oxide_thickness + eventdata->right_oxide_thickness))); // convert voltage to electric field
    if(eventdata->device > 0) eventdata->Set_field(eventdata->voltage); // convert voltage to electric field

    
    std::cout << "graph initialized" << endl;
    std::cout << "max pair degree: " << graph->maxpairdegree() << endl;
    std::cout << "hopping distance: " << graph->hopdist() << endl;
    std::cout << "simulation box size: " << graph->simboxsize() << endl;
    if(eventdata->device > 0) {
        std::cout << "number of left electrode injector nodes " << graph->left()->links().size() << endl;
        std::cout << "number of right electrode injector nodes " << graph->right()->links().size() << endl;
    }
    std::cout << "number of nodes " << graph->Numberofnodes() << endl;    
    
    if(eventdata->device >0) {    
        longrange = new Longrange(graph,eventdata);
        if(eventdata->longrange_slab) longrange->Initialize_slab(graph,eventdata);
        else longrange->Initialize(eventdata);
    }
    else {
        longrange = new Longrange();
    }
    
    std::cout << "longrange profile initialized" << endl;

    state = new StateDevice();
    state->InitStateDevice();
    
    site_inject_probs = new Bsumtree();
    site_inject_probs->initialize(graph->Numberofnodes()); // take care of electrode nodes
    
    //Random charge distribution (assume homogeneous distributed over device)
    int nrcharges;
//    if(eventdata->device>0) nrcharges = (eventdata->voltage)/(2*Pi*eventdata->coulomb_strength*eventdata->simboxsize.x()*eventdata->simboxsize.x())*graph->Numberofnodes();
    if(eventdata->device==0) nrcharges = eventdata->ho_density*graph->Numberofnodes();
    if(eventdata->traj_store) {
        int nrelectrons = eventdata->nr_electrons;
        int nrholes = eventdata->nr_holes;
        if(nrholes != 0) {
            nrcharges = nrholes;
            state->Random_init_injection((int) Hole, nrcharges, site_inject_probs, graph, eventdata, RandomVariable);
        }
        else if(nrelectrons != 0) {
            nrcharges = nrelectrons;
            state->Random_init_injection((int) Electron, nrcharges, site_inject_probs, graph, eventdata, RandomVariable);
        }
    }
    else {
//        state->Random_init_injection((int) Hole, nrcharges, site_inject_probs, graph, eventdata, RandomVariable);
    }
//    std::cout << "initial nrcharges: " << nrcharges << endl;
    
    std::cout << "charges injected" << endl;
    
    if(state->ReservoirEmpty()) state->Grow(eventdata->growsize, eventdata->maxpairdegree);
    
    std::cout << "state initialized" << endl;
    
    non_injection_rates = new Bsumtree();
    left_injection_rates = new Bsumtree();
    right_injection_rates = new Bsumtree();

    std::cout << "binary tree structures initialized" << endl;

    events = new Events();
    events->Init_non_injection_meshes(eventdata);
    events->Initialize_eventvector(graph,state,longrange,eventdata);
    events->Initialize_rates(non_injection_rates, left_injection_rates, right_injection_rates,eventdata);
    if(eventdata->device>0) events->Init_injection_meshes(state, eventdata);
    events->Initialize_after_charge_placement(graph,state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventdata);

    std::cout << "event vectors and meshes initialized" << endl;

    vssmgroup = new Vssmgroup();
    left_chargegroup = new Vssmgroup();
    right_chargegroup = new Vssmgroup();

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
    delete left_chargegroup;
    delete right_chargegroup;
    delete non_injection_rates;
    delete left_injection_rates;
    delete right_injection_rates;
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

    int numcharges_distrib;
    
    if(eventdata->device >0) {
        for(int it= 0; it< eventdata->number_layers; it++ ){
            std::cout << longrange->number_of_nodes(it) << " ";
        }
        std::cout << endl;
    }
    
    std::cout << "total link x distance : " << graph->total_link_distance_x() << endl;
    std::cout << "average hole energy : " << eventdata->avholeenergy << endl;
    std::cout << "disorder strength: " << graph->stddev_hole_node_energy() << endl;
    
    ofstream trajstore;
    char trajfile[100];
    strcpy(trajfile, eventdata->traj_filename.c_str());
    trajstore.open(trajfile);
    
    ofstream vizstore;
    vizstore.open("vizstore");
    

    votca::tools::vec traject;

/*    vector<votca::tools::vec> trajectories;
    trajectories.clear();
    if(eventdata->traj_store) {
        for (int i = 0; i<state->GetCarrierSize(); i++) {
            if(state->GetCarrier(i)->inbox()) {
                trajectories.push_back(votca::tools::vec(0.0,0.0,0.0));
            }
        }
    }*/
    
    vector< vector< vector <double> > > viz_mesh;
    vector< vector< vector <int> > > num_mesh;
    vector<double> layer_cur;
    vector<int> layer_num;
    
    int viz_meshnr_x = eventdata->nx; 
    int viz_meshnr_y = eventdata->ny; 
    int viz_meshnr_z = eventdata->nz;
    
    double viz_size_x = eventdata->simboxsize.x()/viz_meshnr_x;
    double viz_size_y = eventdata->simboxsize.y()/viz_meshnr_y;
    double viz_size_z = eventdata->simboxsize.z()/viz_meshnr_z;
    
    if(eventdata->viz_store) { 
        layer_cur.resize(viz_meshnr_x);
        layer_num.resize(viz_meshnr_x);
        for(int i=0;i<viz_meshnr_x;i++){
            layer_cur[i]=0.0;
            layer_num[i]=0;
        }

        viz_mesh.resize(viz_meshnr_x);
        num_mesh.resize(viz_meshnr_x);
        for(int i = 0;i<viz_meshnr_x;i++) {
            viz_mesh[i].resize(viz_meshnr_y);
            num_mesh[i].resize(viz_meshnr_y);
            for(int j = 0;j<viz_meshnr_y;j++) {
                viz_mesh[i][j].resize(viz_meshnr_z);
                num_mesh[i][j].resize(viz_meshnr_z);
                for(int k = 0; k<viz_meshnr_z;k++) {
                    viz_mesh[i][j][k]=0.0;
                    num_mesh[i][j][k]=0;
                }
            }
        }

        for(int it = 0; it != graph->Numberofnodes(); it++) 
        {
            Node* node = graph->GetNode(it);
            votca::tools::vec position = node->position();
            if(node->type() == (int) NormalNode) {
                int mesh_pos_x = floor(position.x()/viz_size_x);
                int mesh_pos_y = floor(position.y()/viz_size_y);
                int mesh_pos_z = floor(position.z()/viz_size_z);
                num_mesh[mesh_pos_x][mesh_pos_y][mesh_pos_z]++;
                layer_num[mesh_pos_x]++;
            }
        }
    }

    sim_time = 0.0;
    for (long it = 0; it < 2*eventdata->nr_equilsteps + eventdata->nr_timesteps; it++) {
//    for (long it = 0; it < 100; it++) {
        
        if(eventdata->device == 2) {
            // make sure the number of carriers on the left equals
            left_chargegroup->Recompute_injection(left_injection_rates);
            numcharges_distrib = floor(left_chargegroup->totprobsum()+0.5);
            if(numcharges_distrib == 0 && numoutput->holes() == 0) { numcharges_distrib = 1;} // at least one charge per time in the device
            
            while(numcharges_distrib > 0) { // need to fill nodes
                Event* chosencharge = left_chargegroup->Choose_injection_event(events, 0, left_injection_rates, RandomVariable);
                numoutput->Update_ohmic(chosencharge);
                events->On_execute(chosencharge, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventdata);
                left_chargegroup->Recompute_injection(left_injection_rates);
                numcharges_distrib = floor(left_chargegroup->totprobsum()+0.5);

            }

            // make sure the number of carriers on the left equals
            right_chargegroup->Recompute_injection(right_injection_rates);
            numcharges_distrib = floor(right_chargegroup->totprobsum()+0.5);

            // preference given to the major injecting electrode
            while(numcharges_distrib > 0) { // need to fill nodes
                Event* chosencharge = right_chargegroup->Choose_injection_event(events, 1, right_injection_rates, RandomVariable);
                numoutput->Update_ohmic(chosencharge);
                events->On_execute(chosencharge, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventdata);
                right_chargegroup->Recompute_injection(right_injection_rates);
                numcharges_distrib = floor(right_chargegroup->totprobsum()+0.5);

            }
        }

        if(eventdata->device>0) {
            if(ldiv(it, eventdata->steps_update_longrange).rem == 0 && it>0){
                if(eventdata->longrange_slab) longrange->Update_cache_slab(graph,eventdata);
                else                          longrange->Update_cache(eventdata);

    //            double left_oxide_drop = longrange->Calculate_left_oxide_layer(eventdata);
    //            double right_oxide_drop = longrange->Calculate_right_oxide_layer(eventdata);
    //            double in_organic_voltage = eventdata->voltage - left_oxide_drop - right_oxide_drop;
    //            std::cout << "APPLIED VOLTAGE: " << eventdata->voltage << ", VOLTAGE AFTER OXIDE DROP " << in_organic_voltage << endl;

    //            eventdata->Set_field(in_organic_voltage);

                // update extrapolation of potential drop


                events->Recompute_all_events(state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventdata);
            }
        }
        if(eventdata->device == 0) vssmgroup->Recompute_bulk(non_injection_rates);
        if(eventdata->device == 1) vssmgroup->Recompute_device(non_injection_rates, left_injection_rates, right_injection_rates);
        if(eventdata->device == 2) vssmgroup->Recompute_bulk(non_injection_rates);

        double timestep = vssmgroup->Timestep(RandomVariable);
        sim_time += timestep;
        
        Event* chosenevent;
        if(eventdata->device == 0) chosenevent = vssmgroup->Choose_event_bulk(events, non_injection_rates, RandomVariable);
        if(eventdata->device == 1) chosenevent = vssmgroup->Choose_event_device(events, non_injection_rates, left_injection_rates, right_injection_rates, RandomVariable);
        if(eventdata->device == 2) chosenevent = vssmgroup->Choose_event_bulk(events, non_injection_rates, RandomVariable);

        if(it==0) {
            traject = votca::tools::vec(0.0,0.0,0.0);
        }        

        if(eventdata->viz_store) {
            Link* viz_link = chosenevent->link();
            Node* node1 = viz_link->node1();
            Node* node2 = viz_link->node2();
            votca::tools::vec node1_pos = node1->position();
            votca::tools::vec node2_pos = node2->position();

            if(node1->type() == (int) NormalNode && node2->type() == (int) NormalNode) {
                int mesh1_pos_x = floor(node1_pos.x()/viz_size_x);
                int mesh1_pos_y = floor(node1_pos.y()/viz_size_y);
                int mesh1_pos_z = floor(node1_pos.z()/viz_size_z);
                int mesh2_pos_x = floor(node2_pos.x()/viz_size_x);
                int mesh2_pos_y = floor(node2_pos.y()/viz_size_y);
                int mesh2_pos_z = floor(node2_pos.z()/viz_size_z);
                
                if(node1_pos.x() < node2_pos.x()) {

                    if(!(((mesh1_pos_x == mesh2_pos_x)&&(mesh1_pos_y == mesh2_pos_y))&&(mesh1_pos_z == mesh2_pos_z))) {
                        viz_mesh[mesh1_pos_x][mesh1_pos_y][mesh1_pos_z]+=1.0;
                        viz_mesh[mesh2_pos_x][mesh2_pos_y][mesh2_pos_z]+=1.0;
                        layer_cur[mesh1_pos_x] += 1.0;
                        layer_cur[mesh2_pos_x] += 1.0;
                    }
//                    if(layer_cur[mesh_pos_x] <= viz_mesh[mesh_pos_x][mesh_pos_y][mesh_pos_z]) std::cout << viz_mesh[mesh_pos_x][mesh_pos_y][mesh_pos_z] << " " << layer_cur[mesh_pos_x] << endl;
                }
                if(node1_pos.x() > node2_pos.x()) {
                    if(!(((mesh1_pos_x == mesh2_pos_x)&&(mesh1_pos_y == mesh2_pos_y))&&(mesh1_pos_z == mesh2_pos_z))) {
                        viz_mesh[mesh1_pos_x][mesh1_pos_y][mesh1_pos_z]-=1.0;
                        viz_mesh[mesh2_pos_x][mesh2_pos_y][mesh2_pos_z]-=1.0;
                        layer_cur[mesh1_pos_x] -= 1.0;
                        layer_cur[mesh2_pos_x] -= 1.0;
                    }
 //                   if(layer_cur[mesh_pos_x] <= viz_mesh[mesh_pos_x][mesh_pos_y][mesh_pos_z]) std::cout << viz_mesh[mesh_pos_x][mesh_pos_y][mesh_pos_z] << " " << layer_cur[mesh_pos_x] << endl;
                }
            }
            if(it==888888) {
                for(int ix = 0; ix < viz_meshnr_x; ix++) {
                    for(int iy=0; iy < viz_meshnr_y; iy++) {
                        for(int iz=0; iz < viz_meshnr_z; iz++) {
                           double value;
                           double layer_value;
                           if(num_mesh[ix][iy][iz]==0) num_mesh[ix][iy][iz] = 1;
                           value  = viz_mesh[ix][iy][iz]/(1.0*num_mesh[ix][iy][iz]);
                           layer_value = layer_cur[ix]/(1.0*layer_num[ix]);
                           double current = value/layer_value;
                           vizstore << ix << " " << iy << " " << iz << " " << current << endl;
                           std::cout << viz_mesh[ix][iy][iz] << " " << layer_cur[ix] << " " << num_mesh[ix][iy][iz] << endl;
                        }
                    }
                } 
                vizstore.flush();
            }            
        }

        if(eventdata->traj_store) {
            votca::tools::vec traj_hop = chosenevent->link()->r12();
            traject = traject+traj_hop;
        }
        
        if(eventdata->traj_store && ldiv(it,eventdata->nr_reportsteps).rem == 0) {
            trajstore << sim_time << "\t" << traject.x() << "\t" << traject.y() << "\t" << traject.z() << endl;
        }
        
        numoutput->Update(chosenevent, sim_time, timestep); 
        events->On_execute(chosenevent, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventdata);
        
        // check for direct repeats        
        int goto_node_id = chosenevent->link()->node2()->id();
        int from_node_id = chosenevent->link()->node1()->id();
        if(goto_node_id == old_from_node_id && from_node_id == old_to_node_id) repeat_counter++;
        old_from_node_id = from_node_id;
        old_to_node_id = goto_node_id;

        
         if(!eventdata->traj_store &&(it == eventdata->nr_equilsteps || it == 2*eventdata->nr_equilsteps)) numoutput->Init_convergence_check(sim_time);
        
        // equilibration
   
        if(!eventdata->traj_store &&(it == eventdata->nr_equilsteps || it == 2*eventdata->nr_equilsteps)) {
            numoutput->Initialize_equilibrate();
            sim_time = 0.0;
        }
        // convergence checking
        
        if(!eventdata->traj_store && (ldiv(it,10000).rem==0 && it> 2*eventdata->nr_equilsteps)) numoutput->Convergence_check(sim_time, eventdata);

        // direct output
        if(!eventdata->traj_store && ldiv(it,10000).rem==0){
            std::cout << it << " " << repeat_counter << " " << 
                         numoutput->iv_conv() << " " << numoutput->iv_count() << " " << 
                         numoutput->reco_conv() << " " << numoutput->reco_count() <<  " " << 
                         sim_time << " " << timestep << " " << vssmgroup->totprobsum() << " "  << vssmgroup->noninjectprobsum() << " "  << vssmgroup->leftinjectprobsum() << " "  << vssmgroup->rightinjectprobsum() << " " ;
            numoutput->Write(sim_time);
            std::cout << endl;
//            std::cout << "position" << endl;
//           
//            for(int i =0; i< eventdata->number_layers; i++) {
//                std::cout << longrange->position(i) << " ";
//            }
//            
//            std::cout << endl;
            if(eventdata->device>0) {
                std::cout << "charge" << endl;
           
                for(int i =0; i< eventdata->number_layers; i++) {
                    std::cout << longrange->Get_cached_density(i, eventdata) << " ";
                }
            
                std::cout << endl;
            }
//            std::cout << "pot" << endl;
//           
//            for(int i =0; i< eventdata->number_layers; i++) {
//                std::cout << longrange->Get_cached_longrange(i) << " ";
//            }            
//            std::cout << endl;
        }
        
        if(eventdata->traj_store && ldiv(it,10000).rem==0){
            std::cout << "step " << it << " of " << 2*eventdata->nr_equilsteps + eventdata->nr_timesteps << " timesteps" << endl;
        }        
        
        // break out of loop
//        if(numoutput->iv_conv() && numoutput->reco_conv()) {break;}
        if(!eventdata->traj_store && numoutput->iv_conv()) {break;}
        
    }

}

}}


#endif	/* __VOTCA_KMC_DIODE_H */
