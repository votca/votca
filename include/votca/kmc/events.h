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

#ifndef __VOTCA_KMC_EVENTS_H_
#define __VOTCA_KMC_EVENTS_H_


#include <votca/kmc/graph.h>
#include <votca/kmc/state.h>
#include <votca/kmc/event.h>
#include <votca/kmc/bsumtree.h>
#include <votca/kmc/longrange.h>
#include <votca/kmc/globaleventinfo.h>

namespace votca { namespace kmc {
  
using namespace std;

enum action{Add, Remove };

class Events {
    
public:
    
    vector<Event*> El_non_injection_events;
    vector<Event*> Ho_non_injection_events;
    vector<Event*> El_injection_events;
    vector<Event*> Ho_injection_events;
    Bsumtree* El_non_injection_rates;
    Bsumtree* Ho_non_injection_rates;
    Bsumtree* El_injection_rates;
    Bsumtree* Ho_injection_rates;
    Longrange* longrange;
    
    int nholes;
    int nelectrons;
    int ncarriers;
    
    void On_execute(Event* event, Graph* graph, State* state, Globaleventinfo* globevent);

    void Recompute_all_injection_events(Graph* graph, Globaleventinfo* globevent);
    void Recompute_all_non_injection_events(Graph* graph, State* state, Globaleventinfo* globevent);
  
    void Initialize_eventvector(Graph* graph, State* state, Globaleventinfo* globevent);
    void Initialize_longrange(Graph* graph, Globaleventinfo* globevent);
    
    bool el_dirty;
    bool ho_dirty;
    
private:
    void Initialize_injection_eventvector(Node* electrode, vector<Event*> eventvector, CarrierType cartype);
    void Grow_non_injection_eventvector(int carrier_grow_size, vector<Carrier*> carriers, vector<Event*> eventvector,int max_pair_degree);

    void Add_remove_carrier(action AR, Carrier* carrier, Graph* graph, Node* action_node, State* state, Globaleventinfo* globevent);
    void Effect_potential_and_non_injection_rates(action AR, Carrier* carrier, Graph* graph, State* state, Globaleventinfo* globevent);
    void Effect_injection_rates(action AR, Graph* graph, Carrier* carrier, double dist_to_electrode, Node* electrode, Globaleventinfo* globevent);    
    
    double Compute_Coulomb_potential(double startx, myvec dif, myvec sim_box_size, Globaleventinfo* globevent);
};

void Events::Initialize_longrange(Graph* graph, Globaleventinfo* globevent) {
    longrange = new Longrange();
    longrange->Initialize(graph,globevent); 
}

void Events::On_execute(Event* event, Graph* graph, State* state, Globaleventinfo* globevent) {
    
    if(event->fromtype == Fromtransfer) {
        Carrier* carrier = event->carrier;
        Node* fromnode = graph->nodes[carrier->carrier_node_ID]; 
        Node* tonode = fromnode->pairing_nodes[event->tonode_ID];
        Add_remove_carrier(Remove,carrier,graph,fromnode,state,globevent);
    
        if(event->totype == Totransfer) {
            Add_remove_carrier(Add,carrier,graph,tonode,state,globevent);
        }
        else if(event->totype == Recombination) {
            Carrier* recombined_carrier = tonode->carriers_on_node[0];
            Add_remove_carrier(Remove, recombined_carrier, graph, tonode,state,globevent);
            if(carrier->carrier_type == Electron) {
                state->Sell(state->electrons, state->electron_reservoir, carrier->carrier_ID);
                state->Sell(state->holes, state->hole_reservoir, recombined_carrier->carrier_ID);
            }
            else if(carrier->carrier_type == Hole) {
                state->Sell(state->holes, state->hole_reservoir, carrier->carrier_ID);
                state->Sell(state->electrons, state->electron_reservoir, recombined_carrier->carrier_ID);
            }
        }
        else if(event->totype == Collection) {
            if(carrier->carrier_type == Electron) {
                state->Sell(state->electrons, state->electron_reservoir, carrier->carrier_ID);
            }
            else if(carrier->carrier_type == Hole) {
                state->Sell(state->holes, state->hole_reservoir, carrier->carrier_ID);
            }            
        }
    }
    else if(event->fromtype == Injection) {
        Node* tonode = event->electrode->pairing_nodes[event->tonode_ID];
        
        if(event->totype == Totransfer) {
            int carrier_ID;
            if(event->inject_cartype == Electron) {
                if(state->electron_reservoir.empty()){
                    state->Grow(state->electrons, state->electron_reservoir, globevent->state_grow_size, graph->max_pair_degree);
                    Grow_non_injection_eventvector(globevent->state_grow_size, state->electrons, El_non_injection_events,graph->max_pair_degree);
                    El_non_injection_rates->resize(El_non_injection_events.size());
                }
                carrier_ID = state->Buy(state->electrons, state->electron_reservoir);
                Add_remove_carrier(Add,state->electrons[carrier_ID],graph,tonode,state,globevent);
            }
            else if (event->inject_cartype == Hole) {
                if(state->hole_reservoir.empty()){
                    state->Grow(state->holes, state->hole_reservoir, globevent->state_grow_size, graph->max_pair_degree);
                    Grow_non_injection_eventvector(globevent->state_grow_size, state->holes, Ho_non_injection_events,graph->max_pair_degree);
                    Ho_non_injection_rates->resize(Ho_non_injection_events.size());
                }  
                carrier_ID = state->Buy(state->holes, state->hole_reservoir);
                Add_remove_carrier(Add,state->holes[carrier_ID],graph,tonode,state,globevent);
            }
        }
        else if(event->totype == Recombination) {
            Carrier* recombined_carrier = tonode->carriers_on_node[0];
            Add_remove_carrier(Remove,recombined_carrier,graph,tonode,state,globevent);
            if(event->inject_cartype == Electron) {
                state->Sell(state->holes, state->hole_reservoir, recombined_carrier->carrier_ID);
            }
            else if(event->inject_cartype == Hole) {
                state->Sell(state->electrons, state->electron_reservoir, recombined_carrier->carrier_ID);
            }                        
        }
    }
}

void Events::Add_remove_carrier(action AR, Carrier* carrier,Graph* graph, Node* action_node, State* state, Globaleventinfo* globevent){

    if(AR == Add) {
        carrier->carrier_node_ID = action_node->node_ID;
        action_node->carriers_on_node.push_back(carrier);
    
        if (carrier->carrier_type == Hole) {
            nholes++;
        }
        else if (carrier->carrier_type == Electron) {
           nelectrons++;
        } 
        ncarriers++;

        state->Add_to_coulomb_mesh(graph, carrier, globevent);
    }
    else if(AR == Remove) {
        action_node->carriers_on_node.pop_back();
        // Remove existing carrier from lattice
        if (carrier->carrier_type == Hole) {
            nholes--;
        }
        else if (carrier->carrier_type == Electron) {
            nelectrons--;
        }
        ncarriers--;
    }
    
    Effect_potential_and_non_injection_rates(AR,carrier,graph,state, globevent);
 
    // check proximity to left electrode
    if(globevent->device){
        double dist_to_left_electrode = action_node->node_position.x();
        if(dist_to_left_electrode<graph->hopdist){
            Effect_injection_rates(AR,graph,carrier,dist_to_left_electrode,graph->left_electrode,globevent);
        }
    
        // check proximity to right electrode
        double dist_to_right_electrode = graph->sim_box_size.x() - action_node->node_position.x();
        if(dist_to_right_electrode<graph->hopdist){
            Effect_injection_rates(AR,graph,carrier,dist_to_right_electrode,graph->right_electrode, globevent);
        }
    }  
  
    if(AR == Remove){
        state->Remove_from_coulomb_mesh(graph, carrier, globevent);
    }
}


void Events::Effect_potential_and_non_injection_rates(action AR, Carrier* carrier, Graph* graph, State* state,
                                                   Globaleventinfo* globevent) {

    int interact_sign;
    
    if(AR == Add) {interact_sign = 1;}
    if(AR == Remove) {interact_sign = -1;}
    if(carrier->carrier_type == Electron) {interact_sign *= -1;}
    if(carrier->carrier_type == Hole) {interact_sign *=1;}
    
    Node* carnode = graph->nodes[carrier->carrier_node_ID];
    
    //calculate the change to the longrange cache
    
    if(globevent->device){ 
        int layer_index = carnode->layer_index;
        longrange->layercharge[layer_index] += interact_sign;
    }
     
    myvec carpos = carnode->node_position;

    // Define cubic boundaries in non-periodic coordinates
    double ix1 = carpos.x()-globevent->coulcut-graph->hopdist; double ix2 = carpos.x()+globevent->coulcut+graph->hopdist;
    double iy1 = carpos.y()-globevent->coulcut-graph->hopdist; double iy2 = carpos.y()+globevent->coulcut+graph->hopdist;
    double iz1 = carpos.z()-globevent->coulcut-graph->hopdist; double iz2 = carpos.z()+globevent->coulcut+graph->hopdist;

    // Break periodicity in x-direction
    if(globevent->device) {
        if (ix1<0.0) ix1 = 0.0;
        if (ix2>=graph->sim_box_size.x()) ix2 = graph->sim_box_size.x();
    }
  
    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    int sx1 = floor(ix1/globevent->coulcut);
    int sx2 = floor(ix2/globevent->coulcut);
    int sy1 = floor(iy1/globevent->coulcut);
    int sy2 = floor(iy2/globevent->coulcut);
    int sz1 = floor(iz1/globevent->coulcut);
    int sz2 = floor(iz2/globevent->coulcut);
  
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
        int r_isz = isz;
        while (r_isz < 0) r_isz += state->meshsizeZ;
        while (r_isz >= state->meshsizeZ) r_isz -= state->meshsizeZ;
        for (int isy=sy1; isy<=sy2; isy++) {
            int r_isy = isy;
            while (r_isy < 0) r_isy += state->meshsizeY;
            while (r_isy >= state->meshsizeY) r_isy -= state->meshsizeY;
            for (int isx=sx1; isx<=sx2; isx++) {
                int r_isx = isx;
                while (r_isx < 0) r_isx += state->meshsizeX;
                while (r_isx >= state->meshsizeX) r_isx -= state->meshsizeX;
                for (int icartype = 0;icartype <2;icartype++) {
                
                    // Ask a list of all charges in this sublattice
                    list<int>::iterator li1,li2,li3;
                    list<int> *carrierList = &state->coulomb_mesh[r_isx][r_isy][r_isz][icartype];
                    li1 = carrierList->begin();
                    li2 = carrierList->end();
                    for (li3=li1; li3!=li2; li3++) {
                        int probecarrier_ID = *li3;
                        Carrier* probecarrier;
                        if(icartype == 0) {probecarrier = state->electrons[probecarrier_ID];}
                        if(icartype == 1) {probecarrier = state->holes[probecarrier_ID];}
                        Node* probenode = graph->nodes[probecarrier->carrier_node_ID];
                        myvec probepos = probenode->node_position;
                        int probecharge;
                        if(icartype == 0) {
                            probecharge = -1;
                        }
                        else {
                            probecharge = 1;
                        }
          
                        interact_sign *= probecharge;
          
                        // Compute coordinates in non-periodic lattice
                        myvec periodic_convert = myvec((isx-r_isx)*globevent->coulcut,(isy-r_isy)*globevent->coulcut,(isz-r_isz)*globevent->coulcut);
                        myvec np_probepos = probepos + periodic_convert;
                        myvec distance = np_probepos-carpos;

                        double distancesqr = abs(distance)*abs(distance);

                        if (probecarrier_ID!=carrier->carrier_ID) {
                            if((carnode->node_ID!=probenode->node_ID)&&(distancesqr<=globevent->coulcut*globevent->coulcut)) { 
                                
                                // Charge interacting with its own images, taken care off in graph.h
                                // In case multiple charges are on the same node, coulomb calculation on the same spot is catched
                          
                                //First we take the direction sr interactions into account
                                if (AR==Add) carrier->srfrom +=interact_sign*Compute_Coulomb_potential(np_probepos.x(),-1.0*distance,
                                                            graph->sim_box_size,globevent);
                                probecarrier->srfrom += interact_sign*Compute_Coulomb_potential(carpos.x(),distance,
                                                            graph->sim_box_size,globevent);
                            }
                            if (AR==Add) {
              
                                // Adjust Coulomb potential for neighbours of the added carrier
                                for (unsigned int jump=0; jump < carnode->pairing_nodes.size(); jump++) {
                                    myvec jumpdistancevector = carnode->static_event_info[jump].distance;
                                    myvec jumpcarrierpos = carnode->node_position + jumpdistancevector;
                                    myvec jumpdistance = np_probepos - jumpcarrierpos;
                                    double distancejumpsqr = abs(jumpdistance)*abs(jumpdistance);

                                    if(distancejumpsqr <= globevent->coulcut*globevent->coulcut) {
                                    
                                        carrier->srto[jump] += interact_sign*Compute_Coulomb_potential(np_probepos.x(),jumpdistance,
                                                         graph->sim_box_size, globevent);
                                    }
                                }
                            }
                            else if (AR==Remove) {
                            
                                // Reset Coulomb potential for carrier1 and its neighbours
                                carrier->srfrom = 0.0;
                                for (unsigned int jump=0; jump < carnode->pairing_nodes.size(); jump++) {
                                    carrier->srto[jump] = 0.0;  
                                }                            
                            }
           
                            // Adjust Coulomb potential and event rates for neighbours of carrier2
                            for (unsigned int jump=0; jump < probenode->pairing_nodes.size(); jump++) {
                                myvec jumpdistancevector = probenode->static_event_info[jump].distance;
                                myvec jumpprobepos = np_probepos+jumpdistancevector;
                                myvec jumpdistance = carpos-jumpprobepos;
                                double distsqr = abs(jumpdistance)*abs(jumpdistance);
                                int event_ID = probecarrier->carrier_ID*graph->max_pair_degree+jump;
                                
                                double fromlongrange;
                                double tolongrange;
                                if(globevent->device) {
                                    fromlongrange = longrange->Get_cached_longrange(probenode->layer_index);
                                    if(probenode->pairing_nodes[jump]->node_type == Normal) {
                                        tolongrange = longrange->Get_cached_longrange(probenode->pairing_nodes[jump]->layer_index);
                                    }
                                    else { // collection
                                        tolongrange = 0.0;
                                    }
                                }
                                else {
                                    fromlongrange = 0.0;
                                    tolongrange = 0.0;
                                }
                                
                                if(distsqr <= globevent->coulcut*globevent->coulcut) {
                                    if(probecarrier->carrier_type==Electron) {
                                        probecarrier->srto[jump] += 
                                                        interact_sign*Compute_Coulomb_potential(carpos.x(),jumpdistance,
                                                        graph->sim_box_size, globevent);
                                        
                                        El_non_injection_events[event_ID]->Set_non_injection_event(graph->nodes, probecarrier, jump, fromlongrange, tolongrange, globevent);
                                        El_non_injection_rates->setrate(event_ID, El_non_injection_events[event_ID]->rate);
                                        el_dirty = true;
                                    }
                                    else if(probecarrier->carrier_type==Hole) {
                                        probecarrier->srto[jump] += 
                                                            interact_sign*Compute_Coulomb_potential(carpos.x(),jumpdistance,
                                                            graph->sim_box_size, globevent);
                                        Ho_non_injection_events[event_ID]->Set_non_injection_event(graph->nodes, probecarrier, jump, fromlongrange, tolongrange, globevent);
                                        Ho_non_injection_rates->setrate(event_ID, Ho_non_injection_events[event_ID]->rate);
                                        ho_dirty = true;                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }   
    }  

    // update event rates for carrier 1 , done after all carriers within radius coulcut are checked
    for (unsigned int jump=0; jump < carnode->pairing_nodes.size(); jump++) {
        int event_ID = carrier->carrier_ID*graph->max_pair_degree+jump;
    
        double fromlongrange;
        double tolongrange;
        if(globevent->device) {
            fromlongrange = longrange->Get_cached_longrange(carnode->layer_index);
            if(carnode->pairing_nodes[jump]->node_type == Normal) {
                tolongrange = longrange->Get_cached_longrange(carnode->pairing_nodes[jump]->layer_index);
            }
            else { // collection
                tolongrange = 0.0;
            }
        }
        else {
            fromlongrange = 0.0;
            tolongrange = 0.0;
        }
        
        if(carrier->carrier_type==Electron) {
            if(AR == Add) {
                El_non_injection_events[event_ID]->Set_non_injection_event(graph->nodes, carrier, jump, fromlongrange, tolongrange, globevent);
                El_non_injection_rates->setrate(event_ID, El_non_injection_events[event_ID]->rate);
                el_dirty = true;
            }
            else {
                El_non_injection_events[event_ID]->fromtype = Fromnotinbox;
                El_non_injection_events[event_ID]->totype = Tonotinbox;
                El_non_injection_events[event_ID]->rate = 0.0;
                El_non_injection_rates->setrate(event_ID, 0.0);
                el_dirty = true;
            }
        }
        else if(carrier->carrier_type==Hole) {
            if(AR == Add) {
                Ho_non_injection_events[event_ID]->Set_non_injection_event(graph->nodes, carrier, jump, fromlongrange, tolongrange, globevent);
                Ho_non_injection_rates->setrate(event_ID, Ho_non_injection_events[event_ID]->rate);
                ho_dirty = true;
            }
            else {
                Ho_non_injection_events[event_ID]->fromtype = Fromnotinbox;
                Ho_non_injection_events[event_ID]->totype = Tonotinbox;
                Ho_non_injection_events[event_ID]->rate = 0.0;
                Ho_non_injection_rates->setrate(event_ID, 0.0);
                ho_dirty = true;
            }            
        }
    }
}        
        
void Events::Effect_injection_rates(action AR, Graph* graph, Carrier* carrier, 
                                                   double dist_to_electrode, Node* electrode, 
                                                   Globaleventinfo* globevent) {
                                                   
    int interact_sign;
    int x_mesh;
    
    if(AR == Add) {interact_sign = 1;}
    if(AR == Remove) {interact_sign = -1;}
    if(carrier->carrier_type == Electron) {interact_sign *= -1;}
    if(carrier->carrier_type == Hole) {interact_sign *= 1;}
    if(electrode->node_type == LeftElectrode) {x_mesh = 0;}
    if(electrode->node_type == RightElectrode) {x_mesh = graph->nodemeshsizeX-1;}
    
    Node* carnode = graph->nodes[carrier->carrier_node_ID];
    myvec carpos = carnode->node_position;
  
    double bound = sqrt(double(globevent->coulcut*globevent->coulcut - dist_to_electrode*dist_to_electrode));

    // Define cubic boundaries in non-periodic coordinates
    double iy1 = carpos.y()-bound-graph->hopdist; double iy2 = carpos.y()+bound+graph->hopdist;
    double iz1 = carpos.z()-bound-graph->hopdist; double iz2 = carpos.z()+bound+graph->hopdist;

    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    int sy1 = floor(iy1/graph->hopdist);
    int sy2 = floor(iy2/graph->hopdist);
    int sz1 = floor(iz1/graph->hopdist);
    int sz2 = floor(iz2/graph->hopdist);
  
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
        int r_isz = isz;
        while (r_isz < 0) r_isz += graph->nodemeshsizeZ;
        while (r_isz >= graph->nodemeshsizeZ) r_isz -= graph->nodemeshsizeZ;
        for (int isy=sy1; isy<=sy2; isy++) {
            int r_isy = isy;
            while (r_isy < 0) r_isy += graph->nodemeshsizeY;
            while (r_isy >= graph->nodemeshsizeY) r_isy -= graph->nodemeshsizeY;
       
            // Ask a list of all nodes in this sublattice
            list<Node*>::iterator li1,li2,li3;
            list<Node*> *nodemesh = &graph->node_mesh[x_mesh][r_isy][r_isz];
            li1 = nodemesh->begin();
            li2 = nodemesh->end();
            for (li3=li1; li3!=li2; li3++) {
                Node* probenode = *li3;
                myvec probepos = probenode->node_position;
          
                // Compute coordinates in non-periodic lattice
          
                myvec periodic_convert = myvec(0.0,(isy-r_isy)*graph->hopdist,(isz-r_isz)*graph->hopdist);
                myvec np_probepos = probepos + periodic_convert;
                myvec distance = np_probepos-carpos;

                double distancesqr = abs(distance)*abs(distance);

                if ((probenode->node_ID!=carnode->node_ID)&&(distancesqr <= globevent->coulcut*globevent->coulcut)) { // calculated for holes, multiply interact_sign with -1 for electrons
                    probenode->injection_potential +=interact_sign*Compute_Coulomb_potential(carpos.x(),distance,graph->sim_box_size,globevent);
                    int event_ID;
                    int injector_ID;
                    double tolongrange;
                        
                    if(probenode->node_type == Normal){
                        tolongrange = longrange->Get_cached_longrange(probenode->layer_index);
                    }
                    else { // collection (in this case injection to collection)
                        tolongrange = 0.0;
                    }
                        
                    if(electrode->node_type == LeftElectrode) {
                        injector_ID = probenode->left_injector_ID;
                        event_ID = injector_ID;
                        if(globevent->left_injection[1]){
                            Ho_injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Hole, 0.0, tolongrange, globevent);
                            Ho_injection_rates->setrate(event_ID, Ho_injection_events[event_ID]->rate);
                            ho_dirty = true;
                        }
                        if(globevent->left_injection[0]) {
                            El_injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Electron, 0.0, tolongrange, globevent);
                            El_injection_rates->setrate(event_ID, El_injection_events[event_ID]->rate);
                            el_dirty = true;
                        }
                    }
                    else if(electrode->node_type == RightElectrode) {
                        injector_ID = probenode->right_injector_ID;
                        if(globevent->right_injection[1]) {
                            event_ID = injector_ID;
                            if(globevent->left_injection[1]) event_ID += graph->nr_left_injector_nodes;
                            Ho_injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Hole, 0.0, tolongrange, globevent);
                            Ho_injection_rates->setrate(event_ID, Ho_injection_events[event_ID]->rate);
                            ho_dirty = true;
                        }                               
                        if(globevent->right_injection[0]) {
                            event_ID = injector_ID;
                            if(globevent->left_injection[0]) event_ID += graph->nr_left_injector_nodes;
                            El_injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Electron, 0.0, tolongrange, globevent);
                            El_injection_rates->setrate(event_ID, El_injection_events[event_ID]->rate);
                            el_dirty = true;
                        }
                    }
                }
            }
        }
    }  
}

double Events::Compute_Coulomb_potential(double startx, myvec dif, myvec sim_box_size, Globaleventinfo* globevent) {

    double coulpot;
    double RC = globevent->coulcut;
    double RCSQR = RC*RC;
   
    if(!globevent->device) {
        coulpot = 1.0/abs(dif)-1.0/RC;
    }
    else {
        coulpot = 1.0/abs(dif)-1.0/RC; // self image potential is taken into account elsewhere
        
        double L = sim_box_size.x();
        double distsqr_planar = dif.y()*dif.y() + dif.z()*dif.z();
        //double distsqr = dif.x()*dif.x() + distsqr_planar;
      
        int sign;
        double distx_1;
        double distx_2;
        double distancesqr_1;
        double distancesqr_2;
        bool outside_cut_off1 = false;
        bool outside_cut_off2 = false;
      
        while(!(outside_cut_off1&&outside_cut_off2)) {
            for (int i=0;i<globevent->nr_sr_images; i++) {
                if (div(i,2).rem==0) { // even generation
                    sign = -1;
                    distx_1 = i*L + 2*startx + dif.x();
                    distx_2 = (i+2)*L - 2*startx - dif.x(); 
                }
                else {
                    sign = 1;
                    distx_1 = (i+1)*L + dif.x();
                    distx_2 = (i+1)*L - dif.x();
                }
                distancesqr_1 = distx_1*distx_1 + distsqr_planar;
                if (distancesqr_1<=RCSQR) {
                    coulpot += sign*1.0/sqrt(distancesqr_1)-1.0/(RC);
                }
                else {
                    outside_cut_off1 = true;
                }
                distancesqr_2 = distx_2*distx_2 + distsqr_planar;
                if (distancesqr_2<=RCSQR) {
                    coulpot += sign*1.0/sqrt(distancesqr_2)-1.0/(RC);
                }
                else {
                    outside_cut_off2 = true;
                }
            }
        }
    }
    return coulpot;
}


void Events::Recompute_all_non_injection_events(Graph* graph, State* state, Globaleventinfo* globevent) {
    
    int Event_map;
       
    for (unsigned int electron_ID = 0; electron_ID<state->electrons.size(); electron_ID++) {
        Node* electron_node = graph->nodes[state->electrons[electron_ID]->carrier_node_ID];
        for (unsigned int ipair = 0; ipair < electron_node->pairing_nodes.size();ipair++){
            
            Event_map = electron_ID*graph->max_pair_degree + ipair;
            Carrier* electron = state->electrons[electron_ID];
            
            double lrfrom;
            double lrto;
            
            if(globevent->device && electron->is_in_sim_box){
                lrfrom = longrange->Get_cached_longrange(electron_node->layer_index);
                if(electron_node->pairing_nodes[ipair]->node_type == Normal) {
                    lrto = longrange->Get_cached_longrange(electron_node->pairing_nodes[ipair]->layer_index);
                }
                else { // Collection
                    lrto = 0.0;
                }
            }
            else {
                lrfrom = 0.0;
                lrto = 0.0;
            }
            
            El_non_injection_events[Event_map]->Set_non_injection_event(graph->nodes,electron, ipair, lrfrom,lrto, globevent);
            El_non_injection_rates->setrate(Event_map,El_non_injection_events[Event_map]->rate);
            el_dirty = true;
        }
    }

    for (unsigned int hole_ID = 0; hole_ID<state->holes.size(); hole_ID++) {
        Node* hole_node = graph->nodes[state->holes[hole_ID]->carrier_node_ID];
        for (unsigned int ipair = 0; ipair < hole_node->pairing_nodes.size();ipair++){
            
            Event_map = hole_ID*graph->max_pair_degree + ipair;
            Carrier* hole = state->holes[hole_ID];
            
            double lrfrom;
            double lrto;
            
            if(globevent->device && hole->is_in_sim_box){
                lrfrom = longrange->Get_cached_longrange(hole_node->layer_index);
                if(hole_node->pairing_nodes[ipair]->node_type == Normal) {
                    lrto = longrange->Get_cached_longrange(hole_node->pairing_nodes[ipair]->layer_index);
                }
                else { // Collection
                    lrto = 0.0;
                }
            }
            else {
                lrfrom = 0.0;
                lrto = 0.0;
            }
            
            Ho_non_injection_events[Event_map]->Set_non_injection_event(graph->nodes,hole, ipair, lrfrom ,lrto, globevent);
            Ho_non_injection_rates->setrate(Event_map,Ho_non_injection_events[Event_map]->rate);
            ho_dirty = true;
        }
    }
}

void Events::Recompute_all_injection_events(Graph* graph, Globaleventinfo* globevent) {
    
    int Event_map;

    for (unsigned int inject_node = 0; inject_node<graph->left_electrode->pairing_nodes.size(); inject_node++) {

        Event_map = inject_node;

        double lrto;
        if(graph->left_electrode->pairing_nodes[inject_node]->node_type == Normal) {
            lrto = longrange->Get_cached_longrange(graph->left_electrode->pairing_nodes[inject_node]->layer_index);
        }
        else { // Collection
            lrto = 0.0;
        }
        if(globevent->left_injection[0]){
            El_injection_events[Event_map]->Set_injection_event(graph->left_electrode, inject_node, Electron, 0.0, lrto, globevent);   
            El_injection_rates->setrate(Event_map,El_injection_events[Event_map]->rate);
            el_dirty = true;
        }
        if(globevent->left_injection[1]) {
            Ho_injection_events[Event_map]->Set_injection_event(graph->left_electrode, inject_node, Hole, 0.0, lrto, globevent);   
            Ho_injection_rates->setrate(Event_map,Ho_injection_events[Event_map]->rate);
            ho_dirty = true;
        }        
    }
    
    for (unsigned int inject_node = 0; inject_node<graph->right_electrode->pairing_nodes.size(); inject_node++) {

        double lrto;
        if(graph->right_electrode->pairing_nodes[inject_node]->node_type == Normal) {
            lrto = longrange->Get_cached_longrange(graph->right_electrode->pairing_nodes[inject_node]->layer_index);
        }
        else { // Collection
            lrto = 0.0;
        }        
        
        
        if(globevent->right_injection[0]){
            Event_map = inject_node;    
            if(globevent->left_injection[0]) Event_map = graph->nr_left_injector_nodes + inject_node;
            El_injection_events[Event_map]->Set_injection_event(graph->right_electrode, inject_node, Electron, 0.0 , lrto, globevent);
            El_injection_rates->setrate(Event_map,El_injection_events[Event_map]->rate);
            el_dirty = true;
        }
        if(globevent->right_injection[1]){
            Event_map = inject_node;    
            if(globevent->left_injection[0]) Event_map = graph->nr_left_injector_nodes + inject_node;
            Ho_injection_events[Event_map]->Set_injection_event(graph->right_electrode, inject_node, Hole, 0.0, lrto, globevent);
            Ho_injection_rates->setrate(Event_map,Ho_injection_events[Event_map]->rate);
            ho_dirty = true;
        }
    }
}

void Events::Initialize_eventvector(Graph* graph, State* state, Globaleventinfo* globevent){ //
    
    El_non_injection_events.clear();
    Ho_non_injection_events.clear();
    Grow_non_injection_eventvector(state->electrons.size(), state->electrons, El_non_injection_events, graph->max_pair_degree);
    Grow_non_injection_eventvector(state->holes.size(), state->holes,Ho_non_injection_events, graph->max_pair_degree);
    El_non_injection_rates->initialize(El_non_injection_events.size());
    Ho_non_injection_rates->initialize(Ho_non_injection_events.size());
    
    if(globevent->device){
        El_injection_events.clear();
        Ho_injection_events.clear();    
        if(globevent->left_injection[0]) Initialize_injection_eventvector(graph->left_electrode,El_injection_events, Electron);
        if(globevent->left_injection[1]) Initialize_injection_eventvector(graph->left_electrode,Ho_injection_events, Hole);
        if(globevent->right_injection[0]) Initialize_injection_eventvector(graph->right_electrode,El_injection_events, Electron);
        if(globevent->right_injection[1]) Initialize_injection_eventvector(graph->right_electrode,Ho_injection_events, Hole);
        El_injection_rates->initialize(El_injection_events.size());
        Ho_injection_rates->initialize(Ho_injection_events.size());
    }
}

void Events::Initialize_injection_eventvector(Node* electrode, vector<Event*> eventvector, CarrierType cartype){

    for (unsigned int inject_node = 0; inject_node<electrode->pairing_nodes.size(); inject_node++) {

        Event *newEvent = new Event();
        eventvector.push_back(newEvent);
        newEvent->electrode = electrode;
        newEvent->inject_cartype = cartype;
        newEvent->tonode_ID = inject_node;
        
    } 
}

void Events::Grow_non_injection_eventvector(int carrier_grow_size, vector<Carrier*> carriers, vector<Event*> eventvector,int max_pair_degree){
    
    int old_nr_carriers = div(eventvector.size(),max_pair_degree).quot; //what was the number of carriers that we started with?
    
    for(int carrier_ID = old_nr_carriers; carrier_ID<old_nr_carriers+carrier_grow_size; carrier_ID++) {
        for(int jump_ID = 0; jump_ID<max_pair_degree;jump_ID++) {

            Event *newEvent = new Event();
            eventvector.push_back(newEvent);

            newEvent->carrier = carriers[carrier_ID];
            newEvent->tonode_ID = jump_ID;
        }         
    }    
}


}} 

#endif
