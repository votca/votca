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
    }
     
    ~Events() {
        typename std::vector<Event*>::iterator it;
        for (it = _non_injection_events.begin(); it != _non_injection_events.end(); it++ ) delete *it;
        for (it = _injection_events.begin(); it != _injection_events.end(); it++ ) delete *it;
    } 
    
    /// On execute methode
    void On_execute(Event* event, GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo);
    /// On execute method node wise
    void On_execute_node(Node* node, int action, int carrier_type, GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo);
    /// Execute method in case carrier is added on node
    void Add_carrier(Node* node, int carrier_type, GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo);
    /// Execute method in case carrier is removed from node
    void Remove_carrier(Node* node, GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo);
    
    /// Recalculate rates of all events
    void Recompute_all_events(StateDevice* state, Longrange* longrange,Bsumtree* non_injection_rates, Bsumtree* injection_rates,  Eventinfo* eventinfo);
    /// Recalculate rates of all non-injection events
    void Recompute_all_non_injection_events(StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Eventinfo* eventinfo);
    /// Recalculate rates of all injection events
    void Recompute_all_injection_events(StateDevice* state, Longrange* longrange, Bsumtree* injection_rates, Eventinfo* eventinfo);    

    /// Initialize event vectors
    void Initialize_eventvector(GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    /// Initialize injection event vector
    void Initialize_injection_eventvector(int Event_counter, Node* electrode, int carrier_type, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    /// Grow (and initialize) non-injection event vector
    void Grow_non_injection_eventvector(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    
    /// Initialize rates (after initialization of events)
    void Initialize_rates(Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo);
    
    /// Initialize mesh/potential and rates after placement of charges
    void Initialize_after_charge_placement(GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo);
    
    /// Initialize mesh for non-injection events
    void Init_non_injection_meshes(StateDevice* state, Eventinfo* eventinfo);
    /// Initialize mesh for injection events
    void Init_injection_meshes(StateDevice* state, Eventinfo* eventinfo);
    /// Resize mesh
    vector< vector< vector <list<int> > > > Resize_mesh(int meshnr_x, int meshnr_y, int meshnr_z);
    /// Add id (of node or carrier) to mesh
    void Add_to_mesh(int ID, votca::tools::vec position, Eventinfo* eventinfo);
    /// Remove id (of node or carrier) from mesh
    void Remove_from_mesh(int ID,votca::tools::vec position, Eventinfo* eventinfo);

    /// Effect of adding/removing carrier to/from box on the coulomb potential and event rates
    void Effect_potential_and_rates(int action, CarrierDevice* carrier, Node* node, GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo);    
    /// Effect of adding/removing carrier to/from box on the coulomb potential and event rates (non-injection)
    void Effect_potential_and_non_injection_rates(int action, CarrierDevice* carrier1, Node* node, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Eventinfo* eventinfo);
    /// Effect of adding/removing carrier to/from box on the coulomb potential and event rates (injection)
    void Effect_injection_rates(int action, CarrierDevice* carrier, Node* node, Node* electrode, double dist_to_electrode, StateDevice* state, Longrange* longrange, Bsumtree* injection_rates,  Eventinfo* eventinfo);
    
    /// Calculate shortrange coulomb potential
    /// startx is x-coordinate of carrier of which we want to calculate the coulomb potential
    /// dif is startx - x coordinate of node on which we want to know the coulomb potential
    double Compute_Coulomb_potential(double startx, votca::tools::vec dif, bool direct, votca::tools::vec sim_box_size, Eventinfo* eventinfo);

    /// Obtain non-injection event with id eventID
    Event* get_non_injection_event(int eventID) {return _non_injection_events[eventID];}
    /// Obtain injection event with id eventID
    Event* get_injection_event(int eventID) {return _injection_events[eventID];}
    
private:

    vector<Event*> _non_injection_events;
    vector<Event*> _injection_events;
    
    double _meshsize_x; int _inject_meshnr_x; double _meshsize_y; double _meshsize_z; 
    vector< vector< vector <list<int> > > > _non_injection_events_mesh;
    vector< vector< vector <list<int> > > > _left_injection_events_mesh;
    vector< vector< vector <list<int> > > > _right_injection_events_mesh;

};

void Events::On_execute(Event* event, GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo) {
    
    Node* node1 = event->link()->node1();
    Node* node2 = event->link()->node2();
    
    On_execute_node(node1, event->action_node1(), event->carrier_type(), graph, state, longrange, non_injection_rates, injection_rates, eventinfo );
    On_execute_node(node2, event->action_node2(), event->carrier_type(), graph, state, longrange, non_injection_rates, injection_rates, eventinfo );
}

void Events::On_execute_node(Node* node, int action, int carrier_type, GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo) {
        
    if(action == (int) None)        {                                                                                                                   }
    else if(action == (int) Add)    { Add_carrier(node, carrier_type, graph, state, longrange, non_injection_rates, injection_rates, eventinfo);} 
    else if(action == (int) Remove) { Remove_carrier(node, graph, state, longrange, non_injection_rates, injection_rates, eventinfo);           }
}

void Events:: Add_carrier(Node* node, int carrier_type, GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo) {

    
    int new_carrier_ID;

    //make sure the carrier_reservoir is not empty
    if(state->ReservoirEmpty()){
        state->Grow(eventinfo->growsize, eventinfo->maxpairdegree);
        Grow_non_injection_eventvector(state, longrange, eventinfo);
        non_injection_rates->resize(_non_injection_events.size());
    }

    //"buy" the "new" carrier
    new_carrier_ID = state->Buy();
    CarrierDevice* new_carrier = state->GetCarrier(new_carrier_ID);
    
    if(carrier_type == (int) Electron) {
        new_carrier->SetCarrierType((int) Electron);
    }
    else if(carrier_type == (int) Hole) {
        new_carrier->SetCarrierType((int) Hole);
    }

    //place the new carrier in the graph
    new_carrier->SetCarrierNode(node);
    node->AddCarrier(new_carrier_ID);

    //add to mesh
    Add_to_mesh(new_carrier_ID,node->position(), eventinfo);

    //Determine effect on Coulomb potential and rates
    Effect_potential_and_rates((int) Add, new_carrier, node, graph, state, longrange, non_injection_rates, injection_rates, eventinfo);
}

void Events:: Remove_carrier(Node* node, GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo) {
   
    CarrierDevice* removed_carrier = state->GetCarrier(node->occ());

    //remove from graph
    node->RemoveCarrier();    

    //Determine effect on Coulomb potential and rates
    Effect_potential_and_rates((int) Remove, removed_carrier, node, graph, state, longrange, non_injection_rates, injection_rates, eventinfo);    

    //remove from mesh
    Remove_from_mesh(removed_carrier->id(), node->position(), eventinfo);

    //push back to reservoir
    state->Sell(removed_carrier->id());    
}

inline void Events::Effect_potential_and_rates(int action, CarrierDevice* carrier, Node* node, GraphDevice* graph, StateDevice* state, Longrange* longrange, 
                                   Bsumtree* non_injection_rates, Bsumtree* injection_rates,Eventinfo* eventinfo) {

    Effect_potential_and_non_injection_rates(action, carrier, node, state, longrange, non_injection_rates, eventinfo);
    
    if(eventinfo->device){
        votca::tools::vec nodeposition = node->position();
        if(nodeposition.x()<eventinfo->coulcut) {

            double dist_to_electrode = nodeposition.x();
            Effect_injection_rates(action, carrier, node, graph->left(), dist_to_electrode, state, longrange, injection_rates, eventinfo);
        }
        else if(eventinfo->simboxsize.x()-nodeposition.x() < eventinfo->coulcut) {

            double dist_to_electrode = eventinfo->simboxsize.x()-nodeposition.x();
            Effect_injection_rates(action, carrier, node, graph->right(), dist_to_electrode, state, longrange, injection_rates, eventinfo);
        }
    }
}

void Events::Effect_potential_and_non_injection_rates(int action, CarrierDevice* carrier1, Node* node, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Eventinfo* eventinfo) {

    int interact_sign;
    if(action == (int) Add)    {interact_sign =  1;}
    if(action == (int) Remove) {interact_sign = -1;}

    if(carrier1->type() == (int) Electron) {interact_sign *= -1;}
    if(carrier1->type() == (int) Hole)     {interact_sign *=  1;}
   
    double RCSQR = eventinfo->coulcut*eventinfo->coulcut;
    
    if(eventinfo->device){ 
        longrange->Add_charge(interact_sign, dynamic_cast<NodeDevice*>(node)->layer());
    }

    votca::tools::vec carrier1_pos = node->position();

    // Define cubic boundaries in non-periodic coordinates
    double ix1 = carrier1_pos.x()-eventinfo->coulcut-eventinfo->hopdist; double ix2 = carrier1_pos.x()+eventinfo->coulcut+eventinfo->hopdist;
    double iy1 = carrier1_pos.y()-eventinfo->coulcut-eventinfo->hopdist; double iy2 = carrier1_pos.y()+eventinfo->coulcut+eventinfo->hopdist;
    double iz1 = carrier1_pos.z()-eventinfo->coulcut-eventinfo->hopdist; double iz2 = carrier1_pos.z()+eventinfo->coulcut+eventinfo->hopdist;

    // Break periodicity in x-direction
    if(eventinfo->device) {
        if (ix1<0.0) ix1 = 0.0;
        if (ix2>eventinfo->simboxsize.x()) ix2 = eventinfo->simboxsize.x();
    }

    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    int sx1 = floor(ix1/_meshsize_x);
    int sx2 = floor(ix2/_meshsize_x);
    int sy1 = floor(iy1/_meshsize_y);
    int sy2 = floor(iy2/_meshsize_y);
    int sz1 = floor(iz1/_meshsize_z);
    int sz2 = floor(iz2/_meshsize_z);
    
    // Catch for the exceptional case that a carrier is on the boundary of the box
    if(sx2 == eventinfo->mesh_x) {sx2--;}
    if(sy2 == eventinfo->mesh_y) {sy2--;}
    if(sz2 == eventinfo->mesh_z) {sz2--;}
    
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
        int r_isz = isz;
        while (r_isz < 0) r_isz += eventinfo->mesh_z;
        while (r_isz >= eventinfo->mesh_z) r_isz -= eventinfo->mesh_z;
        for (int isy=sy1; isy<=sy2; isy++) {
            int r_isy = isy;
            while (r_isy < 0) r_isy += eventinfo->mesh_y;
            while (r_isy >= eventinfo->mesh_y) r_isy -= eventinfo->mesh_y;
            for (int isx=sx1; isx<=sx2; isx++) {
                int r_isx = isx;
                while (r_isx < 0) r_isx += eventinfo->mesh_x;
                while (r_isx >= eventinfo->mesh_x) r_isx -= eventinfo->mesh_x;
                
                // Ask a list of all charges in this sublattice
                list<int>::iterator li1,li2,li3;
                list<int> *carrierlist = &_non_injection_events_mesh[r_isx][r_isy][r_isz];
                li1 = carrierlist->begin();
                li2 = carrierlist->end();
                for (li3=li1; li3!=li2; li3++) {

                    int carrier2_ID = *li3;
                    CarrierDevice* carrier2 = state->GetCarrier(carrier2_ID);

                    int carrier2_type = carrier2->type();
                    Node* carrier2_node = carrier2->node();
                    votca::tools::vec carrier2_pos = carrier2_node->position();

                    int carrier2_charge;
                    if(carrier2_type == (int) Electron) {
                        carrier2_charge = -1;
                    }
                    else if(carrier2_type == (int) Hole) {
                        carrier2_charge = 1;
                    }

                    interact_sign *= carrier2_charge;
          
                    // Compute coordinates in non-periodic lattice
                    votca::tools::vec periodic_convert = votca::tools::vec((isx-r_isx)*_meshsize_x,(isy-r_isy)*_meshsize_y,(isz-r_isz)*_meshsize_z);
                    votca::tools::vec np_carrier2_pos = carrier2_pos + periodic_convert;
                    votca::tools::vec distance = np_carrier2_pos-carrier1_pos;

                    double distancesqr = abs(distance)*abs(distance);

                    if (carrier2_ID==carrier1->id()) { // self_image_potential interaction taken into account elsewhere (graph.h)
                    }
                    else {
                        if((node->id()!=carrier2_node->id())&&(distancesqr<=RCSQR)) { 

                            // Charge interacting with its own images, taken care off in graph.h
                          
                            //First we take the direct sr interactions into account
                            if (action == (int) Add) carrier1->Add_from_Coulomb(interact_sign*Compute_Coulomb_potential(np_carrier2_pos.x(),-1.0*distance, true, eventinfo->simboxsize,eventinfo));
                            carrier2->Add_from_Coulomb(interact_sign*Compute_Coulomb_potential(carrier1_pos.x(),distance, true, eventinfo->simboxsize,eventinfo));
                        }
                        if (action== (int) Add) {

                            // Adjust Coulomb potential for neighbours of the added carrier
                            for (int it = 0 ; it < node->links().size(); it++) {
                                
                                votca::tools::vec jumpdistancevector = node->links()[it]->r12();
                                votca::tools::vec jump_from_carrier1_pos = carrier1_pos + jumpdistancevector;
                                votca::tools::vec jumpdistance = np_carrier2_pos-jump_from_carrier1_pos;
                                double distancejumpsqr = abs(jumpdistance)*abs(jumpdistance);

                                if(distancejumpsqr <= RCSQR) {
                                    if (node->links()[it]->node2()->id() != carrier2_node->id()) { // coulomb interaction on equal nodes lead to zero distance vector (fixed by exciton binding energy)
                                        carrier1->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(np_carrier2_pos.x(),-1.0*jumpdistance, true, eventinfo->simboxsize, eventinfo),it);
                                    }
                                    else {
                                        carrier1->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(np_carrier2_pos.x(),-1.0*jumpdistance, false, eventinfo->simboxsize, eventinfo),it);
                                    }
                                }
                            }
                        }
                        // Adjust Coulomb potential and event rates for neighbours of carrier2
                        typename std::vector<Link*>::iterator it;
                        for (int it = 0; it < carrier2_node->links().size(); it++) {

                            votca::tools::vec jumpdistancevector = carrier2_node->links()[it]->r12();
                            votca::tools::vec jump_from_carrier2_pos = np_carrier2_pos+jumpdistancevector;
                            votca::tools::vec jumpdistance = jump_from_carrier2_pos - carrier1_pos;
                            double distsqr = abs(jumpdistance)*abs(jumpdistance);
                            if(distsqr <= RCSQR) {

                                int event_ID = carrier2->id()*eventinfo->maxpairdegree + it;
                                if(carrier2_node->links()[it]->node2()->id() == node->id()) {
                                    carrier2->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(carrier1_pos.x(),jumpdistance,false,eventinfo->simboxsize, eventinfo), it);
                                    _non_injection_events[event_ID]->Set_event(carrier2_node->links()[it], carrier2_type, state, longrange, eventinfo); // if carrier2 is actually linking with carrier1
                                }
                                else {
                                    carrier2->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(carrier1_pos.x(),jumpdistance,true,eventinfo->simboxsize, eventinfo), it);
                                    _non_injection_events[event_ID]->Determine_rate(state, longrange, eventinfo);
                                }
                                non_injection_rates->setrate(event_ID,_non_injection_events[event_ID]->rate());
                            }
                        }
                    }
                    if (action == (int) Remove) { // will be done in outside method

                        // Reset Coulomb potential for carrier1 and its neighbours
                        carrier1->Set_from_Coulomb(0.0);
                        carrier1->Reset_to_Coulomb();
                    }
                }
            }
        }   
    }  

    // update event rates for carrier 1 , done after all carriers within radius coulcut are checked
    for (int it = 0; it < node->links().size(); it++) { 
        int event_ID = carrier1->id()*eventinfo->maxpairdegree+it;
    
        if(action == (int) Add) {
            _non_injection_events[event_ID]->Set_event(node->links()[it], carrier1->type(), state, longrange, eventinfo);
            non_injection_rates->setrate(event_ID, _non_injection_events[event_ID]->rate());
        }
        else {
            _non_injection_events[event_ID]->Set_not_in_box_event();
            non_injection_rates->setrate(event_ID, 0.0);
        }
    }
}

void Events:: Effect_injection_rates(int action, CarrierDevice* carrier, Node* node, Node* electrode, double dist_to_electrode, StateDevice* state, Longrange* longrange, Bsumtree* injection_rates, Eventinfo* eventinfo){

    int interact_sign;
    int x_mesh;
    
    if(action == (int) Add) {interact_sign = 1;}
    if(action == (int) Remove) {interact_sign = -1;}
    if(carrier->type() == (int) Electron) {interact_sign *= -1;}
    if(carrier->type() == (int) Hole) {interact_sign *= 1;}

    votca::tools::vec carrier1_pos = node->position();
    double RCSQR = eventinfo->coulcut*eventinfo->coulcut;

    // typically one wants hopdist to be smaller than RC
    double bound;
    if (carrier1_pos.x()>eventinfo->hopdist)       bound = sqrt(double(RCSQR - (dist_to_electrode-eventinfo->hopdist)*(dist_to_electrode-eventinfo->hopdist)));
    else if (carrier1_pos.x()<=eventinfo->hopdist) bound = eventinfo->coulcut;

    // Define cubic boundaries in non-periodic coordinates
    double ix1 = 0.0;                                       double ix2 = eventinfo->hopdist;
    double iy1 = carrier1_pos.y()-bound-eventinfo->hopdist; double iy2 = carrier1_pos.y()+bound+eventinfo->hopdist;
    double iz1 = carrier1_pos.z()-bound-eventinfo->hopdist; double iz2 = carrier1_pos.z()+bound+eventinfo->hopdist;

    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    
    int sx1 = floor(ix1/_meshsize_x); int sx2 = floor(ix2/_meshsize_x);
    int sy1 = floor(iy1/_meshsize_y); int sy2 = floor(iy2/_meshsize_y);
    int sz1 = floor(iz1/_meshsize_z); int sz2 = floor(iz2/_meshsize_z);

    if(sy2 == eventinfo->mesh_y) {sy2--;}
    if(sz2 == eventinfo->mesh_z) {sz2--;}    
   
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
        int r_isz = isz;
        while (r_isz < 0) r_isz += eventinfo->mesh_z;
        while (r_isz >= eventinfo->mesh_z) r_isz -= eventinfo->mesh_z;
        for (int isy=sy1; isy<=sy2; isy++) {
            int r_isy = isy;
            while (r_isy < 0) r_isy += eventinfo->mesh_y;
            while (r_isy >= eventinfo->mesh_y) r_isy -= eventinfo->mesh_y;
            for (int isx = sx1; isx<=sx2; isx++) {
                int r_isx = isx;

                // Ask a list of all nodes in this sublattice
                list<int>::iterator li1,li2,li3;
                list<int> *eventmesh ;

                if(electrode->type() == (int) LeftElectrodeNode) {
                    eventmesh = &_left_injection_events_mesh[r_isx][r_isy][r_isz];
                }
                else if(electrode->type() == (int) RightElectrodeNode) {
                    eventmesh = &_right_injection_events_mesh[r_isx][r_isy][r_isz];
                }

                li1 = eventmesh->begin();
                li2 = eventmesh->end();
                for (li3=li1; li3!=li2; li3++) {
                    int eventID = *li3;
                    Event* inject_event = _injection_events[eventID];
                    Node* probenode = inject_event->link()->node2();
                    votca::tools::vec probepos = probenode->position();
         
                    // Compute coordinates in non-periodic lattice
          
                    votca::tools::vec periodic_convert = votca::tools::vec(0.0,(isy-r_isy)*eventinfo->hopdist,(isz-r_isz)*eventinfo->hopdist);
                    votca::tools::vec np_probepos = probepos + periodic_convert;
                    votca::tools::vec distance = np_probepos-carrier1_pos;

                    double distancesqr = abs(distance)*abs(distance);
                    if(distancesqr <= RCSQR) {

                        if(inject_event->carrier_type() == (int) Electron) {interact_sign *= -1;}
                    
                        if (probenode->id()!=node->id()) {
                            inject_event->Add_injection_potential(interact_sign*Compute_Coulomb_potential(carrier1_pos.x(),distance,true,eventinfo->simboxsize,eventinfo));                    
                            _injection_events[eventID]->Determine_rate(state, longrange, eventinfo);
                            injection_rates->setrate(eventID, _injection_events[eventID]->rate());
                        }
                        else {
                            inject_event->Add_injection_potential(interact_sign*Compute_Coulomb_potential(carrier1_pos.x(),distance,false,eventinfo->simboxsize,eventinfo));                    
                            _injection_events[eventID]->Set_event(inject_event->link(), inject_event->carrier_type(),state, longrange, eventinfo);
                            injection_rates->setrate(eventID, _injection_events[eventID]->rate());
                        }
                    }
                }
            }
        }
    }  
}

double Events::Compute_Coulomb_potential(double startx, votca::tools::vec dif, bool direct, votca::tools::vec simboxsize, Eventinfo* eventinfo) {

    double RC = eventinfo->coulcut;
    double RCSQR = RC*RC;
    double coulpot;
    
    if(direct){
        coulpot = 1.0/abs(dif)-1.0/RC;
    }
    else {
        coulpot = 0.0;
    }
    
    if(eventinfo->device) {
        
        double L = simboxsize.x();
        double distsqr_planar = dif.y()*dif.y() + dif.z()*dif.z();
      
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
                    coulpot += sign*(1.0/sqrt(distancesqr_1)-1.0/(RC));
                }
                else {
                    outside_cut_off1 = true;
                }
                distancesqr_2 = distx_2*distx_2 + distsqr_planar;
                if (distancesqr_2<=RCSQR) {
                    coulpot += sign*(1.0/sqrt(distancesqr_2)-1.0/(RC));
                }
                else {
                    outside_cut_off2 = true;
                }
            }
        }
    }
    return coulpot;
}

void Events::Recompute_all_events(StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo) {
  
    Recompute_all_non_injection_events(state,longrange, non_injection_rates, eventinfo);
    if(eventinfo->device) Recompute_all_injection_events(state,longrange,injection_rates, eventinfo);

}

void Events::Recompute_all_non_injection_events(StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Eventinfo* eventinfo) {
    
    typename std::vector<Event*>::iterator it;
    for(it = _non_injection_events.begin(); it != _non_injection_events.end(); it++) {
        if(((*it)->final_type() != (int) Notinbox)) {
            (*it)->Determine_rate(state, longrange, eventinfo);
        }
        else {
            (*it)->Set_rate(0.0);
        }
        non_injection_rates->setrate((*it)->id(),(*it)->rate());
    }
}

void Events::Recompute_all_injection_events(StateDevice* state, Longrange* longrange, Bsumtree* injection_rates, Eventinfo* eventinfo) {
    
    typename std::vector<Event*>::iterator it;
    for (it = _injection_events.begin(); it!=_injection_events.end(); it++){
        (*it)->Determine_rate(state, longrange, eventinfo);
        injection_rates->setrate((*it)->id(),(*it)->rate());

    }
}


void Events::Initialize_eventvector(GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo){ //

    typename std::vector<Event*>::iterator it;

    _non_injection_events.clear();
    Grow_non_injection_eventvector(state, longrange, eventinfo);

    if(eventinfo->device){
        _injection_events.clear();
        int Event_id_count = 0;
        if(eventinfo->left_electron_injection) {
            Initialize_injection_eventvector(Event_id_count,graph->left(), (int) Electron, state, longrange, eventinfo); 
            Event_id_count += graph->left()->links().size();
        }
        if(eventinfo->left_hole_injection) {
            Initialize_injection_eventvector(Event_id_count,graph->left(), (int) Hole, state, longrange, eventinfo); 
            Event_id_count += graph->left()->links().size();
        }
        if(eventinfo->right_electron_injection) {
            Initialize_injection_eventvector(Event_id_count,graph->right(), (int) Electron, state, longrange, eventinfo); 
            Event_id_count += graph->right()->links().size();
        }
        if(eventinfo->right_hole_injection) {
            Initialize_injection_eventvector(Event_id_count,graph->right(), (int) Hole, state, longrange, eventinfo);
        }
    }
    
}

void Events::Initialize_rates(Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo){

    non_injection_rates->initialize(_non_injection_events.size());
    typename std::vector<Event*>::iterator it; 
    for (it = _non_injection_events.begin(); it!=_non_injection_events.end(); it++) {non_injection_rates->setrate((*it)->id(),(*it)->rate());}     

    if(eventinfo->device) {
        injection_rates->initialize(_injection_events.size());
        for (it = _injection_events.begin(); it!=_injection_events.end(); it++) {injection_rates->setrate((*it)->id(),(*it)->rate());}  
    }
    
}

void Events::Initialize_injection_eventvector(int Event_id_count, Node* electrode, int carrier_type, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo){

    int Event_map = Event_id_count;
    
    for (int it = 0; it < electrode->links().size(); it++) {

        Event *newEvent = new Event(Event_map, electrode->links()[it], carrier_type, state, longrange, eventinfo);
        newEvent->Set_injection_potential(0.0);
        _injection_events.push_back(newEvent);
        Event_map++;
    }
}

void Events::Initialize_after_charge_placement(GraphDevice* graph, StateDevice* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* injection_rates, Eventinfo* eventinfo){

    for(int icar = 0; icar< state->GetCarrierSize(); icar++) {
        CarrierDevice* probe_carrier = state->GetCarrier(icar);
        NodeDevice* carrier_node = dynamic_cast<NodeDevice*>(probe_carrier->node());
        
        if(probe_carrier->inbox()) {
            
           Add_to_mesh(icar,carrier_node->position(), eventinfo);
           Effect_potential_and_rates((int) Add, probe_carrier, carrier_node, graph, state, longrange, non_injection_rates, injection_rates, eventinfo);
        }
    }
}

void Events::Grow_non_injection_eventvector(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo){

    int old_nr_carriers = div(_non_injection_events.size(),eventinfo->maxpairdegree).quot; //what was the number of carriers that we started with?
    
    // if the grow size of the carriers vector is bigger than the growsize defined in eventinfo object, correct for this
    int growth;
    if(eventinfo->growsize < (state->GetCarrierSize()-old_nr_carriers)) {
        growth = (state->GetCarrierSize() - old_nr_carriers);
    }
    else {
        growth = eventinfo->growsize;
    }
    
    for(int carrier_ID = old_nr_carriers; carrier_ID<old_nr_carriers+growth; carrier_ID++) {
        CarrierDevice* probecarrier = state->GetCarrier(carrier_ID);
        int Event_map = carrier_ID*eventinfo->maxpairdegree;
        Event *newEvent;

        for(int ifill = 0; ifill<eventinfo->maxpairdegree; ifill++) {
            newEvent = new Event(Event_map, (int) Notinbox); 
            _non_injection_events.push_back(newEvent); 
            Event_map ++; 
        } // carrier not in box
    }
}

void Events::Init_non_injection_meshes(StateDevice* state, Eventinfo* eventinfo) {

    // determine the dimensions of the meshes
  
    votca::tools::vec simboxsize = eventinfo->simboxsize;
    _meshsize_x = simboxsize.x()/eventinfo->mesh_x; 
    _meshsize_y = simboxsize.y()/eventinfo->mesh_y; 
    _meshsize_z = simboxsize.z()/eventinfo->mesh_z;

    _non_injection_events_mesh = Resize_mesh(eventinfo->mesh_x,eventinfo->mesh_y,eventinfo->mesh_z);

    // initialize meshes
    for (int icar = 0; icar < state->GetCarrierSize(); icar++ ) {
        if(state->GetCarrier(icar)->inbox()){
            CarrierDevice* carrier = state->GetCarrier(icar);
            votca::tools::vec position = carrier->node()->position();

            Add_to_mesh(icar,position,eventinfo);
        }
    }
}

void Events::Init_injection_meshes(StateDevice* state, Eventinfo* eventinfo) {

    _inject_meshnr_x = floor(eventinfo->hopdist/_meshsize_x)+1;

    _left_injection_events_mesh = Resize_mesh(_inject_meshnr_x,eventinfo->mesh_y,eventinfo->mesh_z);
    _right_injection_events_mesh = Resize_mesh(_inject_meshnr_x,eventinfo->mesh_y,eventinfo->mesh_z);    
    
    
    typename std::vector<Event*>::iterator it;
    for (it = _injection_events.begin(); it != _injection_events.end(); it++ ) {
        if((*it)->link()->node1()->type()==LeftElectrodeNode) {
            votca::tools::vec position = (*it)->link()->node2()->position();
            
            double posx = position.x(); int iposx = floor(posx/_meshsize_x);
            double posy = position.y(); int iposy = floor(posy/_meshsize_y);
            double posz = position.z(); int iposz = floor(posz/_meshsize_z);

            if(iposy == eventinfo->mesh_y) {iposy--;}
            if(iposz == eventinfo->mesh_z) {iposz--;}           
            _left_injection_events_mesh[iposx][iposy][iposz].push_back((*it)->id());
        }
        
        //for the righthandside electrode, we take the distance from said electrode

        if((*it)->link()->node1()->type()==RightElectrodeNode) {

            votca::tools::vec simboxsize = eventinfo->simboxsize;
            votca::tools::vec unflipped_position = (*it)->link()->node2()->position();
            votca::tools::vec position = votca::tools::vec(simboxsize.x() - unflipped_position.x(),unflipped_position.y(),unflipped_position.z());

            double posx = position.x(); int iposx = floor(posx/_meshsize_x);
            double posy = position.y(); int iposy = floor(posy/_meshsize_y);
            double posz = position.z(); int iposz = floor(posz/_meshsize_z);

            if(iposy == eventinfo->mesh_y) {iposy--;}
            if(iposz == eventinfo->mesh_z) {iposz--;}
            _right_injection_events_mesh[iposx][iposy][iposz].push_back((*it)->id()); 
        }
    }
}

vector< vector< vector <list<int> > > > Events::Resize_mesh(int meshnr_x, int meshnr_y, int meshnr_z) {
    
    vector< vector< vector <list<int> > > > mesh;
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
    return mesh;
}

void Events::Add_to_mesh(int ID, votca::tools::vec position, Eventinfo* eventinfo){
    
    
    double posx = position.x(); int iposx = floor(posx/_meshsize_x);
    double posy = position.y(); int iposy = floor(posy/_meshsize_y);
    double posz = position.z(); int iposz = floor(posz/_meshsize_z);

    if(iposx == eventinfo->mesh_x) {iposx--;}
    if(iposy == eventinfo->mesh_y) {iposy--;}
    if(iposz == eventinfo->mesh_z) {iposz--;}
    
    _non_injection_events_mesh[iposx][iposy][iposz].push_back(ID);     
};

void Events::Remove_from_mesh(int ID, votca::tools::vec position, Eventinfo* eventinfo){
    
    double posx = position.x(); int iposx = floor(posx/_meshsize_x);
    double posy = position.y(); int iposy = floor(posy/_meshsize_y);
    double posz = position.z(); int iposz = floor(posz/_meshsize_z);

    if(iposx == eventinfo->mesh_x) {iposx--;}
    if(iposy == eventinfo->mesh_y) {iposy--;}
    if(iposz == eventinfo->mesh_z) {iposz--;}
         
    _non_injection_events_mesh[iposx][iposy][iposz].remove(ID);
}

}} 

#endif
