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
    Events() {
        _non_injection_rates = new Bsumtree();
        _injection_rates = new Bsumtree();
    }
     
    ~Events() {
        typename std::vector<Event*>::iterator it;
        for (it = _non_injection_events.begin(); it != _non_injection_events.end(); it++ ) delete *it;
        for (it = _injection_events.begin(); it != _injection_events.end(); it++ ) delete *it;
        delete _non_injection_rates;
        delete _injection_rates;
    } 
    
/*public:

    Longrange* longrange;
    

    
    void On_execute(Event* event, GraphLattice* graph, State* state, Globaleventinfo* globevent);

  
    void Initialize_longrange(GraphLattice* graph, Globaleventinfo* globevent);*/

    void On_execute(Event* event, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo);
    void On_execute_node(Node* node, int action, int carrier_type, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo);
    void Add_carrier(Node* node, int carrier_type, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo);
    void Remove_carrier(Node* node, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo){};
    
    void Recompute_all_non_injection_events(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo);
    void Recompute_all_injection_events(Eventinfo* eventinfo);    

    void Initialize_eventvector(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo);
    void Initialize_injection_eventvector(int Event_counter, Node* electrode, int carrier_type, StateDevice* state, Eventinfo* eventinfo);
    void Grow_non_injection_eventvector(StateDevice* state, int max_pair_degree, Eventinfo* eventinfo);
    
    
    
    void Init_meshes(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph,StateDevice* state, Eventinfo* eventinfo);
    inline void Resize_mesh(int meshnr_x, int meshnr_y, int meshnr_z, vector< vector< vector <list<int> > > > mesh);
    inline void Add_to_mesh(int ID, votca::tools::vec position, vector< vector< vector <list<int> > > > mesh, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, Eventinfo* eventinfo);
    inline void Remove_from_mesh(int ID,votca::tools::vec position,vector< vector< vector <list<int> > > > mesh, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, Eventinfo* eventinfo);
    
    void Effect_potential_and_non_injection_rates(int action, CarrierDevice* carrier, Node* node, StateDevice* state, Eventinfo* eventinfo, double hopdist, votca::tools::vec simboxsize){};
    
//    void Add_remove_carrier(int action_flag1, Carrier* carrier, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, DNode* action_node, State* state, Globaleventinfo* globevent);
//    void Effect_potential_and_non_injection_rates(action AR, Carrier* carrier, GraphLattice* graph, State* state, Globaleventinfo* globevent);
//    void Effect_injection_rates(action AR, GraphLattice* graph, Carrier* carrier, double dist_to_electrode, DNode* electrode, Globaleventinfo* globevent);    
    
    double Compute_Coulomb_potential(double startx, myvec dif, myvec sim_box_size, Eventinfo* eventinfo);
    
private:

    vector<Event*> _non_injection_events;
    vector<Event*> _injection_events;
    
    int _meshnr_x; int _inject_meshnr_x; int _meshnr_y; int _meshnr_z; 
    vector< vector< vector <list<int> > > > _non_injection_events_mesh;
    vector< vector< vector <list<int> > > > _left_injection_events_mesh;
    vector< vector< vector <list<int> > > > _right_injection_events_mesh;

    Bsumtree* _non_injection_rates;
    Bsumtree* _injection_rates;

    int _nholes;
    int _nelectrons;
    int _ncarriers;
};

/*void Events::Initialize_longrange(GraphLattice* graph, Globaleventinfo* globevent) {
    longrange = new Longrange();
    longrange->Initialize(graph,globevent); 
}*/

void Events::On_execute(Event* event, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo) {
    
    Node* node1 = event->link()->node1();
    Node* node2 = event->link()->node2();
    On_execute_node(node1, event->action_node1(), event->carrier_type(), graph, state, eventinfo );
    On_execute_node(node2, event->action_node2(), event->carrier_type(), graph, state, eventinfo );

}

void Events::On_execute_node(Node* node, int action, int carrier_type, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo) {
    
    if(action == (int) None)        {                                            }
    else if(action == (int) Add)    {Add_carrier(node, carrier_type, graph, state, eventinfo); } 
    else if(action == (int) Remove) {Remove_carrier(node, graph, state, eventinfo);            }
}

void Events:: Add_carrier(Node* node, int carrier_type, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo) {
    
    int new_carrier_ID;
    
    //make sure the carrier_reservoir is not empty
    if(state->ReservoirEmpty()){
        state->Grow(eventinfo->growsize);
        Grow_non_injection_eventvector(state, graph->maxpairdegree(),eventinfo);
        _non_injection_rates->resize(_non_injection_events.size());
    }

    //"buy" the "new" carrier
    new_carrier_ID = state->Buy();
    CarrierDevice* new_carrier = state->GetCarrier(new_carrier_ID);    
    if(carrier_type == (int) Electron) {
        new_carrier->SetCarrierType((int) Electron);
        _nelectrons++;
    }
    else if(carrier_type == (int) Hole) {
        new_carrier->SetCarrierType((int) Hole);
        _nholes++;
    }
    _ncarriers++;

    //place the new carrier in the graph
    
    new_carrier->SetCarrierNode(node);
    node->AddCarrier(new_carrier_ID);
    
//    Add_non_injection_event_to_mesh(Event* event,vector< vector< vector <list<int> > > > mesh, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, Eventinfo* eventinfo);
    
//    state->Add_to_coulomb_mesh(graph, carrier, globevent);
    
    Effect_potential_and_non_injection_rates((int) Add, new_carrier, node, state, eventinfo, graph->hopdist(), graph->simboxsize());   
    
}    
    
 
/*void Events::Effect_potential_and_non_injection_rates(int action, CarrierDevice* carrier, Node* node, StateDevice* state, Eventinfo* eventinfo, double hopdist, votca::tools::vec simboxsize) {

    int interact_sign;
    
    if(action == (int) Add)    {interact_sign =  1;}
    if(action == (int) Remove) {interact_sign = -1;}
    if(carrier->type() == (int) Electron) {interact_sign *= -1;}
    if(carrier->type() == (int) Hole)     {interact_sign *=  1;}
    
    double RCSQR = eventinfo->coulcut*eventinfo->coulcut;
    
    //calculate the change to the longrange cache
    
//    if(globevent->device){ 
//        int layer_index = carnode->layer_index;
//        longrange->layercharge[layer_index] += interact_sign;
//    }
     
    votca::tools::vec carpos = node->position();

    // Define cubic boundaries in non-periodic coordinates
    double ix1 = carpos.x()-eventinfo->mesh_x-hopdist; double ix2 = carpos.x()+eventinfo->mesh_x+hopdist;
    double iy1 = carpos.y()-eventinfo->mesh_y-hopdist; double iy2 = carpos.y()+eventinfo->mesh_y+hopdist;
    double iz1 = carpos.z()-eventinfo->mesh_z-hopdist; double iz2 = carpos.z()+eventinfo->mesh_z+hopdist;

    // Break periodicity in x-direction
    if(eventinfo->device) {
        if (ix1<0.0) ix1 = 0.0;
        if (ix2>simboxsize.x()) ix2 = simboxsize.x();
    }

    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    int sx1 = floor(ix1/eventinfo->mesh_x);
    int sx2 = floor(ix2/eventinfo->mesh_x);
    int sy1 = floor(iy1/eventinfo->mesh_y);
    int sy2 = floor(iy2/eventinfo->mesh_y);
    int sz1 = floor(iz1/eventinfo->mesh_z);
    int sz2 = floor(iz2/eventinfo->mesh_z);
  
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
        int r_isz = isz;
        while (r_isz < 0) r_isz += _meshnr_z;
        while (r_isz >= _meshnr_z) r_isz -= _meshnr_z;
        for (int isy=sy1; isy<=sy2; isy++) {
            int r_isy = isy;
            while (r_isy < 0) r_isy += _meshnr_y;
            while (r_isy >= _meshnr_y) r_isy -= _meshnr_y;
            for (int isx=sx1; isx<=sx2; isx++) {
                int r_isx = isx;
                while (r_isx < 0) r_isx += _meshnr_x;
                while (r_isx >= _meshnr_x) r_isx -= _meshnr_x;
                
                // Ask a list of all charges in this sublattice
                list<int>::iterator li1,li2,li3;
                list<int> *carrierlist = &_non_injection_events_mesh[r_isx][r_isy][r_isz];
                li1 = carrierlist->begin();
                li2 = carrierlist->end();
                for (li3=li1; li3!=li2; li3++) {
                    int probecarrier_ID = *li3;
                    Carrier* probecarrier = state->GetCarrier(probecarrier_ID);
                    int probecartype = probecarrier->type();
                    Node* probenode = probecarrier->node();
                    myvec probepos = probenode->position();

                    int probecharge;
                    if(probecartype == (int) Electron) {
                        probecharge = -1;
                    }
                    else if(probecartype == (int) Hole) {
                        probecharge = 1;
                    }
          
                    interact_sign *= probecharge;
          
                    // Compute coordinates in non-periodic lattice
                    myvec periodic_convert = myvec((isx-r_isx)*eventinfo->mesh_x,(isy-r_isy)*eventinfo->mesh_y,(isz-r_isz)*eventinfo->mesh_z);
                    myvec np_probepos = probepos + periodic_convert;
                    myvec distance = np_probepos-carpos;

                    double distancesqr = abs(distance)*abs(distance);

                    if (probecarrier_ID==carrier->id()) {
                        
                        
                    }
                    else {
                        if((node->id()!=probenode->id())&&(distancesqr<=RCSQR)) { 
                                
                            // Charge interacting with its own images, taken care off in graph.h
                          
                            //First we take the direct sr interactions into account
                            if (action == (int) Add) carrier->srfrom +=interact_sign*Compute_Coulomb_potential(np_probepos.x(),-1.0*distance, simboxsize,eventinfo);
                            probecarrier->srfrom += interact_sign*Compute_Coulomb_potential(carpos.x(),distance, simboxsize,eventinfo);
                            
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
    
}*/

/* 
void Events::Add_remove_carrier(action AR, CarrierDevice* carrier,GraphLattice* graph, DNode* action_node, State* state, Globaleventinfo* globevent){

    if(AR == Add) {

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


void Events::Effect_potential_and_non_injection_rates(action AR, CarrierDevice* carrier, GraphLattice* graph, State* state,
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
                    CarrierDevice* probecarrier = state->carriers[probecarrier_ID];
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
}*/

double Events::Compute_Coulomb_potential(double startx, myvec dif, myvec sim_box_size, Eventinfo* eventinfo) {

    double RC = eventinfo->coulcut;
    double RCSQR = RC*RC;

    double coulpot = 1.0/abs(dif)-1.0/RC;
    
    if(eventinfo->device) {
        
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
            for (int i=0;i<eventinfo->nr_sr_images; i++) {
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

void Events::Recompute_all_non_injection_events(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo) {
    
    typename std::vector<Event*>::iterator it;
    for(it = _non_injection_events.begin(); it != _non_injection_events.end(); it++) {
        if(((*it)->final_type() != (int) Notinbox)&&((*it)->final_type() != (int) Notingraph)) {
            (*it)->Determine_rate(eventinfo);
        }
        else {
            (*it)->Set_rate(0.0);
        }
        _non_injection_rates->setrate((*it)->id(),(*it)->rate());
    }
}

void Events::Recompute_all_injection_events(Eventinfo* eventinfo) {
    
    typename std::vector<Event*>::iterator it;
    for (it = _injection_events.begin(); it!=_injection_events.end(); it++){
        (*it)->Determine_rate(eventinfo);
        _injection_rates->setrate((*it)->id(),(*it)->rate());

    }
}


void Events::Initialize_eventvector(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo){ //

    typename std::vector<Event*>::iterator it;
    
    _non_injection_events.clear();
    Grow_non_injection_eventvector(state,graph->maxpairdegree(), eventinfo);
    int size = _non_injection_events.size();
    std::cout << size << endl;
    _non_injection_rates->initialize(size);
    std::cout << size << endl;

    
    for (it = _non_injection_events.begin(); it!=_non_injection_events.end(); it++) {_non_injection_rates->setrate((*it)->id(),(*it)->rate());}    
    
    if(eventinfo->device){
        _injection_events.clear();
        int Event_id_count = 0;
        if(eventinfo->left_electron_injection) {
            Initialize_injection_eventvector(Event_id_count,graph->left(), (int) Electron, state, eventinfo); 
            Event_id_count += graph->left()->links().size();
        }
        if(eventinfo->left_hole_injection) {
            Initialize_injection_eventvector(Event_id_count,graph->left(), (int) Hole, state, eventinfo); 
            Event_id_count += graph->left()->links().size();
        }
        if(eventinfo->right_electron_injection) {
            Initialize_injection_eventvector(Event_id_count,graph->right(), (int) Electron, state, eventinfo); 
            Event_id_count += graph->right()->links().size();
        }
        if(eventinfo->right_hole_injection) {
            Initialize_injection_eventvector(Event_id_count,graph->right(), (int) Hole, state, eventinfo);
        }
        _injection_rates->initialize(_injection_events.size());
        for (it = _injection_events.begin(); it!=_injection_events.end(); it++) {_injection_rates->setrate((*it)->id(),(*it)->rate());}  
    }
    
}

void Events::Initialize_injection_eventvector(int Event_id_count, Node* electrode, int carrier_type, StateDevice* state, Eventinfo* eventinfo){

    int Event_map = Event_id_count;
    vector<Link*> injectorlinks = electrode->links();
    typename std::vector<Link*>::iterator it;    
    for (it = injectorlinks.begin(); it != injectorlinks.end(); it++) {
        Event *newEvent = new Event(Event_map, (*it), carrier_type, eventinfo, state);
        _injection_events.push_back(newEvent);
        Event_map++;
        std::cout << Event_map << endl;
    }
}


void Events::Grow_non_injection_eventvector(StateDevice* state, int max_pair_degree, Eventinfo* eventinfo){
    
    int old_nr_carriers = div(_non_injection_events.size(),max_pair_degree).quot; //what was the number of carriers that we started with?
    
    for(int carrier_ID = old_nr_carriers; carrier_ID<old_nr_carriers+eventinfo->growsize; carrier_ID++) {
        CarrierDevice* probecarrier = state->GetCarrier(carrier_ID);
        int Event_map = carrier_ID*max_pair_degree;
        if(probecarrier->inbox()) { 
            int fillcount = 0;
            vector<Link*> carrierlinks = probecarrier->node()->links();
            typename std::vector<Link*>::iterator it;    
            for(it = carrierlinks.begin(); it != carrierlinks.end();it++){ 
                Event_map += (*it)->id(); 
                Event *newEvent = new Event(Event_map, (*it), probecarrier->type(), eventinfo, state); 
                _non_injection_events.push_back(newEvent); fillcount++; 
            }
            for(int ifill = fillcount; ifill<max_pair_degree; ifill++) { 
                Event_map += ifill; 
                Event *newEvent = new Event(Event_map, (int) Notingraph); 
                _non_injection_events.push_back(newEvent); // non-existent event
            }
        }
        else {
            for(int ifill = 0; ifill<max_pair_degree; ifill++) { 
                Event_map += ifill; 
                Event *newEvent = new Event(Event_map, (int) Notinbox); 
                _non_injection_events.push_back(newEvent); } // carrier not in box
        }
    }
}

void Events::Init_meshes(GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, StateDevice* state, Eventinfo* eventinfo) {

    // determine the dimensions of the meshes
    
    votca::tools::vec simboxsize = graph->simboxsize();
    _meshnr_x = ceil(simboxsize.x()/eventinfo->mesh_x); 
    _meshnr_y = ceil(simboxsize.y()/eventinfo->mesh_y); 
    _meshnr_z = ceil(simboxsize.z()/eventinfo->mesh_z);
    
    _inject_meshnr_x = ceil(graph->hopdist()/eventinfo->mesh_x);
    
    // resize the meshes
    
    Resize_mesh(_meshnr_x,_meshnr_y,_meshnr_z,_non_injection_events_mesh);
    Resize_mesh(_inject_meshnr_x,_meshnr_y,_meshnr_z,_left_injection_events_mesh);
    Resize_mesh(_inject_meshnr_x,_meshnr_y,_meshnr_z,_right_injection_events_mesh);
    
    // initialize meshes
    for (int icar = 0; icar < state->GetCarrierSize(); icar++ ) {
        if(!state->GetCarrier(icar)->inbox()){
            CarrierDevice* carrier = state->GetCarrier(icar);
            votca::tools::vec position = carrier->node()->position();
            Add_to_mesh(icar,position,_non_injection_events_mesh, graph, eventinfo);
        }
    }
    
    typename std::vector<Event*>::iterator it;
    for (it = _injection_events.begin(); it != _injection_events.end(); it++ ) {
        if((*it)->link()->node1()->type()==LeftElectrodeNode) Add_to_mesh((*it)->id(),(*it)->link()->node2()->position(),_left_injection_events_mesh, graph, eventinfo);
        
        //for the righthandside electrode, we take the distance from said electrode

        if((*it)->link()->node1()->type()==RightElectrodeNode) {
            votca::tools::vec simboxsize = graph->simboxsize();
            votca::tools::vec eventpos = votca::tools::vec(simboxsize.x(),0.0,0.0)-(*it)->link()->node2()->position();
            Add_to_mesh((*it)->id(),eventpos,_right_injection_events_mesh, graph, eventinfo);
        }
    }
    
}

inline void Events::Resize_mesh(int meshnr_x, int meshnr_y, int meshnr_z, vector< vector< vector <list<int> > > > mesh) {
    mesh.resize(meshnr_x);
    for(int i = 0;i<meshnr_x;i++) {
        mesh[i].resize(meshnr_y);
        for(int j = 0;j<meshnr_y;j++) {
            mesh[i][j].resize(meshnr_z);
            for(int k = 0; k<meshnr_z;k++) {
                mesh[i][j][k].clear();
            }
        }
    }
}

inline void Events::Add_to_mesh(int ID, votca::tools::vec position, vector< vector< vector <list<int> > > > mesh, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, Eventinfo* eventinfo){
    
    double posx = position.x(); int iposx = floor(posx/eventinfo->mesh_x);
    double posy = position.y(); int iposy = floor(posy/eventinfo->mesh_y);
    double posz = position.z(); int iposz = floor(posz/eventinfo->mesh_z);
     
    mesh[iposx][iposy][iposz].push_back(ID);     
};

inline void Events::Remove_from_mesh(int ID, votca::tools::vec position, vector< vector< vector <list<int> > > > mesh, GraphDevice<GraphSQL, NodeSQL, LinkSQL>* graph, Eventinfo* eventinfo){
    
    double posx = position.x(); int iposx = floor(posx/eventinfo->mesh_x);
    double posy = position.y(); int iposy = floor(posy/eventinfo->mesh_y);
    double posz = position.z(); int iposz = floor(posz/eventinfo->mesh_z);
     
    mesh[iposx][iposy][iposz].remove(ID);
}

}} 

#endif
