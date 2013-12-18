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
    
    void On_execute(Event* event, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    void On_execute_node(Node* node, int action, int carrier_type, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    void Add_carrier(Node* node, int carrier_type, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    void Remove_carrier(Node* node, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    
    void Recompute_all_events(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    void Recompute_all_non_injection_events(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    void Recompute_all_injection_events(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);    

    void Initialize_eventvector(GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    void Initialize_injection_eventvector(int Event_counter, Node* electrode, int carrier_type, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    void Grow_non_injection_eventvector(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    
    
    
    void Init_meshes(StateDevice* state, Eventinfo* eventinfo);
    inline void Resize_mesh(int meshnr_x, int meshnr_y, int meshnr_z, vector< vector< vector <list<int> > > > mesh);
    inline void Add_to_mesh(int ID, votca::tools::vec position, vector< vector< vector <list<int> > > > mesh, Eventinfo* eventinfo);
    inline void Remove_from_mesh(int ID,votca::tools::vec position,vector< vector< vector <list<int> > > > mesh, Eventinfo* eventinfo);

    void Effect_potential_and_rates(int action, CarrierDevice* carrier, Node* node, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);    
    void Effect_potential_and_non_injection_rates(int action, CarrierDevice* carrier1, Node* node, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    void Effect_injection_rates(int action, CarrierDevice* carrier, Node* node, Node* electrode, double dist_to_electrode, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    
    double Compute_Coulomb_potential(double startx, votca::tools::vec dif, votca::tools::vec sim_box_size, Eventinfo* eventinfo);
 
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

void Events::Recompute_all_events(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
    Recompute_all_non_injection_events(state,longrange,eventinfo);
    Recompute_all_injection_events(state,longrange,eventinfo);
}

void Events::On_execute(Event* event, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
    
    Node* node1 = event->link()->node1();
    Node* node2 = event->link()->node2();
    On_execute_node(node1, event->action_node1(), event->carrier_type(), graph, state, longrange, eventinfo );
    On_execute_node(node2, event->action_node2(), event->carrier_type(), graph, state, longrange, eventinfo );

}

void Events::On_execute_node(Node* node, int action, int carrier_type, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
    
    if(action == (int) None)        {                                            }
    else if(action == (int) Add)    {Add_carrier(node, carrier_type, graph, state, longrange, eventinfo); } 
    else if(action == (int) Remove) {Remove_carrier(node, graph, state, longrange, eventinfo);            }
}

void Events:: Add_carrier(Node* node, int carrier_type, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
    
    int new_carrier_ID;
    
    //make sure the carrier_reservoir is not empty
    if(state->ReservoirEmpty()){
        state->Grow(eventinfo->growsize, eventinfo->maxpairdegree);
        Grow_non_injection_eventvector(state, longrange, eventinfo);
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
    
    //add to mesh
    Add_to_mesh(new_carrier_ID,node->position(),_non_injection_events_mesh, eventinfo);
    
    Effect_potential_and_rates((int) Add, new_carrier, node, graph, state, longrange, eventinfo);

}

void Events:: Remove_carrier(Node* node, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
 
    CarrierDevice* removed_carrier = state->GetCarrier(node->occ());

    // Remove existing carrier from lattice
    if (removed_carrier->type() == Hole) {
        _nholes--;
    }
    else if (removed_carrier->type() == Electron) {
        _nelectrons--;
    }
    _ncarriers--;

    Effect_potential_and_rates((int) Remove, removed_carrier, node, graph, state, longrange, eventinfo);    
    
    //remove from mesh
    Remove_from_mesh(removed_carrier->id(), node->position(), _non_injection_events_mesh,  eventinfo);

    //push to reservoir
    state->Sell(removed_carrier->id());    
    
    //remove from graph
    node->RemoveCarrier();
}

inline void Events::Effect_potential_and_rates(int action, CarrierDevice* carrier, Node* node, GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {

    Effect_potential_and_non_injection_rates((int) Add, carrier, node, state, longrange, eventinfo);
    
    votca::tools::vec nodeposition = node->position();
    if(nodeposition.x()<eventinfo->coulcut) {
        double dist_to_electrode = nodeposition.x();
        Effect_injection_rates((int) Add, carrier, node, graph->left(), dist_to_electrode, state, longrange, eventinfo);
    }
    else if(eventinfo->simboxsize.x()-nodeposition.x() < eventinfo->coulcut) {
        double dist_to_electrode = eventinfo->simboxsize.x()-nodeposition.x();
        Effect_injection_rates((int) Add, carrier, node, graph->right(), dist_to_electrode, state, longrange, eventinfo);
    }
}

void Events::Effect_potential_and_non_injection_rates(int action, CarrierDevice* carrier1, Node* node, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {

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
    double ix1 = carrier1_pos.x()-eventinfo->mesh_x-eventinfo->hopdist; double ix2 = carrier1_pos.x()+eventinfo->mesh_x+eventinfo->hopdist;
    double iy1 = carrier1_pos.y()-eventinfo->mesh_y-eventinfo->hopdist; double iy2 = carrier1_pos.y()+eventinfo->mesh_y+eventinfo->hopdist;
    double iz1 = carrier1_pos.z()-eventinfo->mesh_z-eventinfo->hopdist; double iz2 = carrier1_pos.z()+eventinfo->mesh_z+eventinfo->hopdist;

    // Break periodicity in x-direction
    if(eventinfo->device) {
        if (ix1<0.0) ix1 = 0.0;
        if (ix2>eventinfo->simboxsize.x()) ix2 = eventinfo->simboxsize.x();
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
                    votca::tools::vec periodic_convert = votca::tools::vec((isx-r_isx)*eventinfo->mesh_x,(isy-r_isy)*eventinfo->mesh_y,(isz-r_isz)*eventinfo->mesh_z);
                    votca::tools::vec np_carrier2_pos = carrier2_pos + periodic_convert;
                    votca::tools::vec distance = np_carrier2_pos-carrier1_pos;

                    double distancesqr = abs(distance)*abs(distance);

                    if (carrier2_ID==carrier1->id()) { // self_image_potential interaction
                    }
                    else {
                        if((node->id()!=carrier2_node->id())&&(distancesqr<=RCSQR)) { 
                                
                            // Charge interacting with its own images, taken care off in graph.h
                          
                            //First we take the direct sr interactions into account
                            if (action == (int) Add) carrier1->Add_from_Coulomb(interact_sign*Compute_Coulomb_potential(np_carrier2_pos.x(),-1.0*distance, eventinfo->simboxsize,eventinfo));
                            carrier2->Add_from_Coulomb(interact_sign*Compute_Coulomb_potential(carrier1_pos.x(),distance, eventinfo->simboxsize,eventinfo));
                            
                        }
                        if (action== (int) Add) {
              
                            // Adjust Coulomb potential for neighbours of the added carrier
                            for (int it = 0 ; it < node->links().size(); it++) {
                                votca::tools::vec jumpdistancevector = node->links()[it]->r12();
                                votca::tools::vec jump_from_carrier1_pos = carrier1_pos + jumpdistancevector;
                                votca::tools::vec jumpdistance = np_carrier2_pos-jump_from_carrier1_pos;
                                double distancejumpsqr = abs(jumpdistance)*abs(jumpdistance);

                                if(distancejumpsqr <= RCSQR) {
                                    carrier1->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(np_carrier2_pos.x(),jumpdistance, eventinfo->simboxsize, eventinfo),it);
                                }                                
                                
                            }
                        }
                        else if (action == (int) Remove) { // will be done in outside method
                            // Reset Coulomb potential for carrier1 and its neighbours
                            carrier1->Set_from_Coulomb(0.0);
                            carrier1->Reset_to_Coulomb();
                        }
           
                        // Adjust Coulomb potential and event rates for neighbours of carrier2
                        typename std::vector<Link*>::iterator it;
                        for (int it = 0; it < carrier2_node->links().size(); it++) {
                            votca::tools::vec jumpdistancevector = node->links()[it]->r12();
                            votca::tools::vec jump_from_carrier2_pos = np_carrier2_pos+jumpdistancevector;
                            votca::tools::vec jumpdistance = carrier1_pos-jump_from_carrier2_pos;
                            double distsqr = abs(jumpdistance)*abs(jumpdistance);
                            
                            if(distsqr <= RCSQR) {
                                carrier2->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(carrier1_pos.x(),jumpdistance,eventinfo->simboxsize, eventinfo), it);
                                int event_ID = carrier2->id()*eventinfo->maxpairdegree + it;
                                if(carrier2_node->links()[it]->node2()->id() == node->id()) {
                                    _non_injection_events[event_ID]->Set_event(carrier2_node->links()[it], carrier2_type, state, longrange, eventinfo); // if carrier2 is actually linking with carrier1
                                }
                                else {
                                    _non_injection_events[event_ID]->Determine_rate(state, longrange, eventinfo);
                                }
                                _non_injection_rates->setrate(event_ID,_non_injection_events[event_ID]->rate());
                            }
                        }
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
            _non_injection_rates->setrate(event_ID, _non_injection_events[event_ID]->rate());
        }
        else {
            _non_injection_events[event_ID]->Set_not_in_box_event();
            _non_injection_rates->setrate(event_ID, 0.0);
        }
    }
}

void Events:: Effect_injection_rates(int action, CarrierDevice* carrier, Node* node, Node* electrode, double dist_to_electrode, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo){

                                                   
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
    if (eventinfo->coulcut>eventinfo->hopdist)       bound = sqrt(double(RCSQR - (dist_to_electrode-eventinfo->hopdist)*(dist_to_electrode-eventinfo->hopdist)));
    else if (eventinfo->coulcut<=eventinfo->hopdist) bound = eventinfo->coulcut;

    // Define cubic boundaries in non-periodic coordinates
    double ix1 = 0.0;                                       double ix2 = eventinfo->hopdist;
    double iy1 = carrier1_pos.y()-bound-eventinfo->hopdist; double iy2 = carrier1_pos.y()+bound+eventinfo->hopdist;
    double iz1 = carrier1_pos.z()-bound-eventinfo->hopdist; double iz2 = carrier1_pos.z()+bound+eventinfo->hopdist;

    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    int sx1 = floor(ix1/eventinfo->mesh_x); int sx2 = floor(ix2/eventinfo->mesh_x);
    int sy1 = floor(iy1/eventinfo->mesh_y); int sy2 = floor(iy2/eventinfo->mesh_y);
    int sz1 = floor(iz1/eventinfo->mesh_z); int sz2 = floor(iz2/eventinfo->mesh_z);
  
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
        int r_isz = isz;
        while (r_isz < 0) r_isz += _meshnr_z;
        while (r_isz >= _meshnr_z) r_isz -= _meshnr_z;
        for (int isy=sy1; isy<=sy2; isy++) {
            int r_isy = isy;
            while (r_isy < 0) r_isy += _meshnr_y;
            while (r_isy >= _meshnr_y) r_isy -= _meshnr_y;
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

                    if ((probenode->id()!=node->id())&&(distancesqr <= RCSQR)) { // calculated for holes, multiply interact_sign with -1 for electrons
                        if(inject_event->carrier_type() == (int) Electron) interact_sign *= -1;
                        inject_event->Add_injection_potential(interact_sign*Compute_Coulomb_potential(carrier1_pos.x(),distance,eventinfo->simboxsize,eventinfo));
                        _injection_events[eventID]->Determine_rate(state, longrange, eventinfo);
                        _injection_rates->setrate(eventID, _injection_events[eventID]->rate());
                    }
                    else {
                        _injection_events[eventID]->Set_event(inject_event->link(), inject_event->carrier_type(),state, longrange, eventinfo);
                        _injection_rates->setrate(eventID, _injection_events[eventID]->rate());
                    }
                }
            }
        }
    }  
}

double Events::Compute_Coulomb_potential(double startx, votca::tools::vec dif, votca::tools::vec simboxsize, Eventinfo* eventinfo) {

    double RC = eventinfo->coulcut;
    double RCSQR = RC*RC;

    double coulpot = 1.0/abs(dif)-1.0/RC;
    
    if(eventinfo->device) {
        
        double L = simboxsize.x();
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

void Events::Recompute_all_non_injection_events(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
    
    typename std::vector<Event*>::iterator it;
    for(it = _non_injection_events.begin(); it != _non_injection_events.end(); it++) {
        if(((*it)->final_type() != (int) Notinbox)&&((*it)->final_type() != (int) Notingraph)) {
            (*it)->Determine_rate(state, longrange, eventinfo);
        }
        else {
            (*it)->Set_rate(0.0);
        }
        _non_injection_rates->setrate((*it)->id(),(*it)->rate());
    }
}

void Events::Recompute_all_injection_events(StateDevice* state, Longrange* longrange,  Eventinfo* eventinfo) {
    
    typename std::vector<Event*>::iterator it;
    for (it = _injection_events.begin(); it!=_injection_events.end(); it++){
        (*it)->Determine_rate(state, longrange, eventinfo);
        _injection_rates->setrate((*it)->id(),(*it)->rate());

    }
}


void Events::Initialize_eventvector(GraphDevice* graph, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo){ //

    typename std::vector<Event*>::iterator it;
   
    _non_injection_events.clear();
    Grow_non_injection_eventvector(state, longrange, eventinfo);
    _non_injection_rates->initialize(_non_injection_events.size());

    for (it = _non_injection_events.begin(); it!=_non_injection_events.end(); it++) {_non_injection_rates->setrate((*it)->id(),(*it)->rate());}    

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
        _injection_rates->initialize(_injection_events.size());
        for (it = _injection_events.begin(); it!=_injection_events.end(); it++) {_injection_rates->setrate((*it)->id(),(*it)->rate());}  
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


void Events::Grow_non_injection_eventvector(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo){

    int old_nr_carriers = div(_non_injection_events.size(),eventinfo->maxpairdegree).quot; //what was the number of carriers that we started with?
    
    for(int carrier_ID = old_nr_carriers; carrier_ID<old_nr_carriers+eventinfo->growsize; carrier_ID++) {
        CarrierDevice* probecarrier = state->GetCarrier(carrier_ID);
        int Event_map = carrier_ID*eventinfo->maxpairdegree;
        if(probecarrier->inbox()) { 
            int fillcount = 0;
            for(int it = 0; it < probecarrier->node()->links().size();it++){ 
                Link* probelink = probecarrier->node()->links()[it];
                Event_map += probelink->id(); 
                Event *newEvent = new Event(Event_map, probelink, probecarrier->type(), state, longrange, eventinfo); 
                _non_injection_events.push_back(newEvent); fillcount++; 
            }
            for(int ifill = fillcount; ifill<eventinfo->maxpairdegree; ifill++) { 
                Event_map += ifill; 
                Event *newEvent = new Event(Event_map, (int) Notingraph); 
                _non_injection_events.push_back(newEvent); // non-existent event
            }
        }
        else {
            for(int ifill = 0; ifill<eventinfo->maxpairdegree; ifill++) { 
                Event_map += ifill; 
                Event *newEvent = new Event(Event_map, (int) Notinbox); 
                _non_injection_events.push_back(newEvent); } // carrier not in box
        }
    }
}

void Events::Init_meshes(StateDevice* state, Eventinfo* eventinfo) {

    // determine the dimensions of the meshes
    
    votca::tools::vec simboxsize = eventinfo->simboxsize;
    _meshnr_x = ceil(simboxsize.x()/eventinfo->mesh_x); 
    _meshnr_y = ceil(simboxsize.y()/eventinfo->mesh_y); 
    _meshnr_z = ceil(simboxsize.z()/eventinfo->mesh_z);
    
    _inject_meshnr_x = ceil(eventinfo->hopdist/eventinfo->mesh_x);
    
    // resize the meshes
    
    Resize_mesh(_meshnr_x,_meshnr_y,_meshnr_z,_non_injection_events_mesh);
    Resize_mesh(_inject_meshnr_x,_meshnr_y,_meshnr_z,_left_injection_events_mesh);
    Resize_mesh(_inject_meshnr_x,_meshnr_y,_meshnr_z,_right_injection_events_mesh);
    
    // initialize meshes
    for (int icar = 0; icar < state->GetCarrierSize(); icar++ ) {
        if(!state->GetCarrier(icar)->inbox()){
            CarrierDevice* carrier = state->GetCarrier(icar);
            votca::tools::vec position = carrier->node()->position();
            Add_to_mesh(icar,position,_non_injection_events_mesh, eventinfo);
        }
    }
    
    typename std::vector<Event*>::iterator it;
    for (it = _injection_events.begin(); it != _injection_events.end(); it++ ) {
        if((*it)->link()->node1()->type()==LeftElectrodeNode) Add_to_mesh((*it)->id(),(*it)->link()->node2()->position(),_left_injection_events_mesh, eventinfo);
        
        //for the righthandside electrode, we take the distance from said electrode

        if((*it)->link()->node1()->type()==RightElectrodeNode) {
            votca::tools::vec simboxsize = eventinfo->simboxsize;
            votca::tools::vec eventpos = votca::tools::vec(simboxsize.x(),0.0,0.0)-(*it)->link()->node2()->position();
            Add_to_mesh((*it)->id(),eventpos,_right_injection_events_mesh, eventinfo);
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

inline void Events::Add_to_mesh(int ID, votca::tools::vec position, vector< vector< vector <list<int> > > > mesh, Eventinfo* eventinfo){
    
    double posx = position.x(); int iposx = floor(posx/eventinfo->mesh_x);
    double posy = position.y(); int iposy = floor(posy/eventinfo->mesh_y);
    double posz = position.z(); int iposz = floor(posz/eventinfo->mesh_z);
     
    mesh[iposx][iposy][iposz].push_back(ID);     
};

inline void Events::Remove_from_mesh(int ID, votca::tools::vec position, vector< vector< vector <list<int> > > > mesh,  Eventinfo* eventinfo){
    
    double posx = position.x(); int iposx = floor(posx/eventinfo->mesh_x);
    double posy = position.y(); int iposy = floor(posy/eventinfo->mesh_y);
    double posz = position.z(); int iposz = floor(posz/eventinfo->mesh_z);
     
    mesh[iposx][iposy][iposz].remove(ID);
}

}} 

#endif
