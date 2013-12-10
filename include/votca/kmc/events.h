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


#include <votca/kmc/graphdevice.h>
#include <votca/kmc/statedevice.h>
#include <votca/kmc/event.h>
#include <votca/kmc/bsumtree.h>
#include <votca/kmc/longrange.h>
#include <votca/kmc/eventinfo.h>

namespace votca { namespace kmc {
  
using namespace std;

class Events {

public: 
    Events() {}
     
    ~Events() {
        typename std::vector<Event*>::iterator it;
        for (it = _non_injection_events.begin(); it != _non_injection_events.end(); it++ ) delete *it;
        for (it = _injection_events.begin(); it != _injection_events.end(); it++ ) delete *it;
    } 
    
/*public:
    Bsumtree* Non_injection_rates;
    Bsumtree* Injection_rates;
    Longrange* longrange;
    
    int nholes;
    int nelectrons;
    int ncarriers;
    
    void On_execute(Event* event, GraphLattice* graph, State* state, Globaleventinfo* globevent);

  
    void Initialize_longrange(GraphLattice* graph, Globaleventinfo* globevent);*/

    void Recompute_all_non_injection_events(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo);
    void Recompute_all_injection_events(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo);    

    void Initialize_eventvector(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo);
    void Initialize_injection_eventvector(Node* electrode, int carrier_type, StateDevice* state, Eventinfo* eventinfo);
    void Grow_non_injection_eventvector(StateDevice* state, int max_pair_degree, Eventinfo* eventinfo);
    
    int meshsizeX; int meshsizeY; int meshsizeZ;
    void Init_non_injection_mesh(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, Eventinfo* eventinfo);
    void Init_injection_mesh(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, Eventinfo* eventinfo);
    void Add_to_non_injection_mesh(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, Event* event, Eventinfo* eventinfo);
    void Remove_from_non_injection_mesh(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, Event* event, Eventinfo* eventinfo);
   
//    void Add_remove_carrier(int action_flag1, Carrier* carrier, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, DNode* action_node, State* state, Globaleventinfo* globevent);
//    void Effect_potential_and_non_injection_rates(action AR, Carrier* carrier, GraphLattice* graph, State* state, Globaleventinfo* globevent);
//    void Effect_injection_rates(action AR, GraphLattice* graph, Carrier* carrier, double dist_to_electrode, DNode* electrode, Globaleventinfo* globevent);    
    
    double Compute_Coulomb_potential(double startx, myvec dif, myvec sim_box_size, Eventinfo* eventinfo);
    
private:

    vector<Event*> _non_injection_events;
    vector<Event*> _injection_events;
    
    vector< vector< vector <list<int> > > > _non_injection_event_mesh;
    vector< vector <list<int> > >  _injection_event_mesh;
};

/*void Events::Initialize_longrange(GraphLattice* graph, Globaleventinfo* globevent) {
    longrange = new Longrange();
    longrange->Initialize(graph,globevent); 
}

void Events::On_execute(Event* event, GraphLattice* graph, State* state, Globaleventinfo* globevent) {
    
    if(event->fromtype == Fromtransfer) {
        Carrier* carrier = event->carrier;
        DNode* fromnode = graph->nodes[carrier->carrier_node_ID]; 
        DNode* tonode = fromnode->pairing_nodes[event->tonode_ID];
        Add_remove_carrier(Remove,carrier,graph,fromnode,state,globevent);
    
        if(event->totype == Totransfer) {
            Add_remove_carrier(Add,carrier,graph,tonode,state,globevent);
        }
        else if(event->totype == Recombination) {
            Carrier* recombined_carrier = tonode->carriers_on_node[0];
            Add_remove_carrier(Remove, recombined_carrier, graph, tonode,state,globevent);
            state->Sell(carrier->carrier_ID);
            state->Sell(recombined_carrier->carrier_ID);
        }
        else if(event->totype == Collection) {
            state->Sell(carrier->carrier_ID);
        }
    }
    else if(event->fromtype == Injection) {
        DNode* tonode = event->electrode->pairing_nodes[event->tonode_ID];
        
        if(event->totype == Totransfer) {
            int carrier_ID;
            if(state->carrier_reservoir.empty()){
                state->Grow(globevent->state_grow_size, graph->max_pair_degree);
                Grow_non_injection_eventvector(globevent->state_grow_size, state->carriers,graph->max_pair_degree);
                Non_injection_rates->resize(Non_injection_events.size());
            }

            carrier_ID = state->Buy();
            if(event->inject_cartype == Electron) {
                state->carriers[carrier_ID]->carrier_type == Electron;
            }
            else if(event->inject_cartype == Hole) {
                state->carriers[carrier_ID]->carrier_type == Hole;
            }
            
            Add_remove_carrier(Add,state->carriers[carrier_ID],graph,tonode,state,globevent);
        }
        else if(event->totype == Recombination) {
            Carrier* recombined_carrier = tonode->carriers_on_node[0];
            Add_remove_carrier(Remove,recombined_carrier,graph,tonode,state,globevent);
            state->Sell(recombined_carrier->carrier_ID);
        }
    }
}

void Events::Add_remove_carrier(action AR, Carrier* carrier,GraphLattice* graph, DNode* action_node, State* state, Globaleventinfo* globevent){

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


void Events::Effect_potential_and_non_injection_rates(action AR, Carrier* carrier, GraphLattice* graph, State* state,
                                                   Globaleventinfo* globevent) {

    int interact_sign;
    
    if(AR == Add) {interact_sign = 1;}
    if(AR == Remove) {interact_sign = -1;}
    if(carrier->carrier_type == Electron) {interact_sign *= -1;}
    if(carrier->carrier_type == Hole) {interact_sign *=1;}
    
    DNode* carnode = graph->nodes[carrier->carrier_node_ID];
    
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
                
                // Ask a list of all charges in this sublattice
                list<int>::iterator li1,li2,li3;
                list<int> *carrierList = &state->coulomb_mesh[r_isx][r_isy][r_isz];
                li1 = carrierList->begin();
                li2 = carrierList->end();
                for (li3=li1; li3!=li2; li3++) {
                    int probecarrier_ID = *li3;
                    Carrier* probecarrier = state->carriers[probecarrier_ID];
                    CarrierType probecartype = probecarrier->carrier_type;
                    DNode* probenode = graph->nodes[probecarrier->carrier_node_ID];
                    myvec probepos = probenode->node_position;

                    int probecharge;
                    if(probecartype == Electron) {
                        probecharge = -1;
                    }
                    else if(probecartype == Hole) {
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
                            for (int jump=0; jump < carnode->pairing_nodes.size(); jump++) {
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
                            for (int jump=0; jump < carnode->pairing_nodes.size(); jump++) {
                                carrier->srto[jump] = 0.0;  
                            }                            
                        }
           
                        // Adjust Coulomb potential and event rates for neighbours of carrier2
                        for (int jump=0; jump < probenode->pairing_nodes.size(); jump++) {
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
                                probecarrier->srto[jump] += 
                                                   interact_sign*Compute_Coulomb_potential(carpos.x(),jumpdistance,
                                                   graph->sim_box_size, globevent);
                                        
                                Non_injection_events[event_ID]->Set_non_injection_event(graph->nodes, probecarrier, jump, fromlongrange, tolongrange, globevent);
                                Non_injection_rates->setrate(event_ID, Non_injection_events[event_ID]->rate);
                            }
                        }
                    }
                }
            }
        }   
    }  

    // update event rates for carrier 1 , done after all carriers within radius coulcut are checked
    for (int jump=0; jump < carnode->pairing_nodes.size(); jump++) {
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
        
        if(AR == Add) {
            Non_injection_events[event_ID]->Set_non_injection_event(graph->nodes, carrier, jump, fromlongrange, tolongrange, globevent);
            Non_injection_rates->setrate(event_ID, Non_injection_events[event_ID]->rate);
        }
        else {
            Non_injection_events[event_ID]->fromtype = Fromnotinbox;
            Non_injection_events[event_ID]->totype = Tonotinbox;
            Non_injection_events[event_ID]->rate = 0.0;
            Non_injection_rates->setrate(event_ID, 0.0);
        }
    }
}        
        
void Events::Effect_injection_rates(action AR, GraphLattice* graph, Carrier* carrier, 
                                                   double dist_to_electrode, DNode* electrode, 
                                                   Globaleventinfo* globevent) {
                                                   
    int interact_sign;
    int x_mesh;
    
    if(AR == Add) {interact_sign = 1;}
    if(AR == Remove) {interact_sign = -1;}
    if(carrier->carrier_type == Electron) {interact_sign *= -1;}
    if(carrier->carrier_type == Hole) {interact_sign *= 1;}
    if(electrode->node_type == LeftElectrode) {x_mesh = 0;}
    if(electrode->node_type == RightElectrode) {x_mesh = graph->nodemeshsizeX-1;}
    
    DNode* carnode = graph->nodes[carrier->carrier_node_ID];
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
            list<DNode*>::iterator li1,li2,li3;
            list<DNode*> *nodemesh = &graph->node_mesh[x_mesh][r_isy][r_isz];
            li1 = nodemesh->begin();
            li2 = nodemesh->end();
            for (li3=li1; li3!=li2; li3++) {
                DNode* probenode = *li3;
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
                    
                    int left_injection_count = 0;
                    int right_injection_count = 0;
                    if(electrode->node_type == LeftElectrode) {
                        injector_ID = probenode->left_injector_ID;
                        event_ID = injector_ID;
                        if(globevent->left_injection[0]){
                            event_ID += left_injection_count*graph->nr_left_injector_nodes + right_injection_count*graph->nr_right_injector_nodes;
                            
                            Injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Electron, 0.0, tolongrange, globevent);
                            Injection_rates->setrate(event_ID, Injection_events[event_ID]->rate);
                            left_injection_count++;
                        }
                        if(globevent->left_injection[1]) {
                            event_ID += left_injection_count*graph->nr_left_injector_nodes + right_injection_count*graph->nr_right_injector_nodes;
                            
                            Injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Hole, 0.0, tolongrange, globevent);
                            Injection_rates->setrate(event_ID, Injection_events[event_ID]->rate);
                            left_injection_count++;
                        }
 
                    }
                    else if(electrode->node_type == RightElectrode) {
                        injector_ID = probenode->right_injector_ID;
                        event_ID = injector_ID;
                        if(globevent->right_injection[0]) {
                            event_ID += left_injection_count*graph->nr_left_injector_nodes + right_injection_count*graph->nr_right_injector_nodes;

                            Injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Electron, 0.0, tolongrange, globevent);
                            Injection_rates->setrate(event_ID, Injection_events[event_ID]->rate);
                            right_injection_count++;
                        }
                        if(globevent->right_injection[1]) {
                            event_ID += left_injection_count*graph->nr_left_injector_nodes + right_injection_count*graph->nr_right_injector_nodes;
                            
                            Injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Hole, 0.0, tolongrange, globevent);
                            Injection_rates->setrate(event_ID, Injection_events[event_ID]->rate);
                            right_injection_count++;
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
        double distsqr = dif.x()*dif.x() + distsqr_planar;
      
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

*/

void Events::Recompute_all_non_injection_events(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo) {
    
    int Event_map;
       
    for (int carrier_ID = 0; carrier_ID<state->GetCarrierSize() ; carrier_ID++) {
        
        Carrier* probecarrier = state->GetCarrier(carrier_ID);
        
        if(probecarrier->inbox()) {
            
            Node* probenode = probecarrier->node();
            vector<Link*> carrierlinks = probenode->links();
        
            int fillcount = 0;
            typename std::vector<Link*>::iterator it;            
            for (it = carrierlinks.begin(); it != carrierlinks.end(); it++){
            
                Event_map = carrier_ID*graph->maxpairdegree() + (*it)->id();
            
//            double lrfrom;
//            double lrto;
            
//            if(eventinfo->device && carrier->is_in_sim_box){
//                lrfrom = longrange->Get_cached_longrange(carrier_node->layer_index);
//                if(carrier_node->pairing_nodes[ipair]->node_type == Normal) {
//                    lrto = longrange->Get_cached_longrange(carrier_node->pairing_nodes[ipair]->layer_index);
//                }
//                else { // Collection
//                    lrto = 0.0;
//                }
//            }
//            else {
//                lrfrom = 0.0;
//                lrto = 0.0;
//            }
            
                _non_injection_events[Event_map]->Set_event((*it), probecarrier->type(), state, eventinfo);
//            Non_injection_rates->setrate(Event_map,Non_injection_events[Event_map]->rate);
                fillcount++;
            }
            for (int ifill = fillcount; ifill<graph->maxpairdegree(); ifill++) { Event_map = carrier_ID*graph->maxpairdegree() + ifill; _non_injection_events[Event_map]->Set_not_in_box_event();}
        }
        else  {
            for (int ifill = 0; ifill<graph->maxpairdegree(); ifill++) { Event_map = carrier_ID*graph->maxpairdegree() + ifill; _non_injection_events[Event_map]->Set_not_in_box_event(); }
        }
    }
}

void Events::Recompute_all_injection_events(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo) {
    
    Node* left_electrode = graph->left();
    Node* right_electrode = graph->right();
    vector<Link*> left_electrode_links = left_electrode->links();
    vector<Link*> right_electrode_links = right_electrode->links();
    
    typename std::vector<Link*>::iterator it;    
    for (it = left_electrode_links.begin(); it != left_electrode_links.end(); it++) {

//        double lrto;
//        if(graph->left_electrode->pairing_nodes[inject_node]->node_type == Normal) {
//            lrto = longrange->Get_cached_longrange(graph->left_electrode->pairing_nodes[inject_node]->layer_index);
//        }
//        else { // Collection
//            lrto = 0.0;
//        }
        int event_ID = (*it)->id();
        
        if(eventinfo->left_electron_injection){
            _injection_events[event_ID]->Set_event((*it), (int) Electron, state , eventinfo);   
//            Injection_rates->setrate(Event_ID,Injection_events[Event_ID]->rate);
        }
        if(eventinfo->left_hole_injection) {
            if(eventinfo->left_electron_injection) event_ID += left_electrode_links.size();
            _injection_events[event_ID]->Set_event((*it), (int) Hole, state, eventinfo);   
//            Injection_rates->setrate(Event_ID,Injection_events[Event_ID]->rate);
        }        
    }
    
    for (it = right_electrode_links.begin(); it != right_electrode_links.end(); it++) {

//        double lrto;
//        if(graph->right_electrode->pairing_nodes[inject_node]->node_type == Normal) {
//            lrto = longrange->Get_cached_longrange(graph->right_electrode->pairing_nodes[inject_node]->layer_index);
//        }
//        else { // Collection
//            lrto = 0.0;
//        }        
 
        int event_ID = (*it)->id();        
        
        if(eventinfo->right_electron_injection){
            if(eventinfo->left_electron_injection) event_ID += left_electrode_links.size();
            if(eventinfo->left_hole_injection) event_ID += left_electrode_links.size();
            _injection_events[event_ID]->Set_event((*it), (int) Electron, state, eventinfo);
//            Injection_rates->setrate(Event_ID,Injection_events[Event_ID]->rate);
        }
        if(eventinfo->right_hole_injection){
            if(eventinfo->left_electron_injection) event_ID += left_electrode_links.size();
            if(eventinfo->left_hole_injection) event_ID += left_electrode_links.size();
            if(eventinfo->right_electron_injection) event_ID += right_electrode_links.size();
            _injection_events[event_ID]->Set_event((*it), (int) Hole, state, eventinfo);
//            Injection_rates->setrate(Event_ID,Injection_events[Event_ID]->rate);
        }
    }
}


void Events::Initialize_eventvector(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo){ //
    
    _non_injection_events.clear();
//    Grow_non_injection_eventvector(state->carriers.size(), state->carriers, graph->max_pair_degree);
//    Non_injection_rates->initialize(Non_injection_events.size());
    
    if(eventinfo->device){
        _injection_events.clear();
        if(eventinfo->left_electron_injection) Initialize_injection_eventvector(graph->left(), (int) Electron, state, eventinfo);
        if(eventinfo->left_hole_injection) Initialize_injection_eventvector(graph->left(), (int) Hole, state, eventinfo);
        if(eventinfo->right_electron_injection) Initialize_injection_eventvector(graph->right(), (int) Electron, state, eventinfo);
        if(eventinfo->right_hole_injection) Initialize_injection_eventvector(graph->right(), (int) Hole, state, eventinfo);
//        Injection_rates->initialize(Injection_events.size());
    }
}

void Events::Initialize_injection_eventvector(Node* electrode, int carrier_type, StateDevice* state, Eventinfo* eventinfo){

    vector<Link*> injectorlinks = electrode->links();
    typename std::vector<Link*>::iterator it;    
    for (it = injectorlinks.begin(); it != injectorlinks.end(); it++) {
        Event *newEvent = new Event((*it), carrier_type, eventinfo, state);
        _injection_events.push_back(newEvent);
    } 
}


void Events::Grow_non_injection_eventvector(StateDevice* state, int max_pair_degree, Eventinfo* eventinfo){
    
    int old_nr_carriers = div(_non_injection_events.size(),max_pair_degree).quot; //what was the number of carriers that we started with?
    
    for(int carrier_ID = old_nr_carriers; carrier_ID<old_nr_carriers+eventinfo->growsize; carrier_ID++) {
        Carrier* probecarrier = state->GetCarrier(carrier_ID);
        if(probecarrier->inbox()) { 
            int fillcount = 0;
            vector<Link*> carrierlinks = probecarrier->node()->links();
            typename std::vector<Link*>::iterator it;    
            for(it = carrierlinks.begin(); it != carrierlinks.end();it++)      { Event *newEvent = new Event((*it), probecarrier->type(), eventinfo, state); _non_injection_events.push_back(newEvent); fillcount++; }
            for(int ifill = fillcount; ifill<max_pair_degree; ifill++) { Event *newEvent = new Event(); _non_injection_events.push_back(newEvent); } // non-existent event
        }
        else {
            for(int ifill = 0; ifill<max_pair_degree; ifill++) { Event *newEvent = new Event(); _non_injection_events.push_back(newEvent); } // carrier not in box
        }
    }    
}


}} 

#endif
