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

#include <votca/xtp/events.h>


using namespace std;

namespace votca {
    namespace xtp {
        
void Events::On_execute(Event* event, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, 
                          Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo) {

    if(event->action_pair() == (int) Transfer) {
        On_execute_pair(event->link(),event->carrier_type(), graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);
    }
    else {    
      Node* node1 = event->link()->node1();
      Node* node2 = event->link()->node2();
      
      On_execute_node(node1, event->action_node1(), event->carrier_type(), graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);
      On_execute_node(node2, event->action_node2(), event->carrier_type(), graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);

    }
}

void Events::On_execute_pair(Link* link, int carrier_type, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, 
                              Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo) {
        
    Node* node1 = link->node1();
    Node* node2 = link->node2();
    //remove from graph
    CarrierBulk* transfer_carrier = state->GetCarrier(node1->occ());

    transfer_carrier->IncDistance(link->r12());

    node1->RemoveCarrier();    

    Effect_potential_and_rates((int) Remove, transfer_carrier, node1, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);    

    if(eventinfo->coulomb_strength>0.0) Remove_from_mesh(transfer_carrier->id(), node1->position(), eventinfo);        

    //place the carrier at the new position in the graph

    transfer_carrier->SetCarrierNode(node2);

    node2->AddCarrier(transfer_carrier->id());

    if(eventinfo->coulomb_strength>0.0) Add_to_mesh(transfer_carrier->id(),node2->position(), eventinfo);

    Effect_potential_and_rates((int) Add, transfer_carrier, node2, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);    

}

void Events::On_execute_node(Node* node, int action, int carrier_type, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, 
                              Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo) {
        
    if(action == (int) None)        {                                                                                                                   }
    else if(action == (int) Add)    { Add_carrier(node, carrier_type, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);} 
    else if(action == (int) Remove) { Remove_carrier(node, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);           }
}

void Events:: Add_carrier(Node* node, int carrier_type, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, 
                            Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo) {

    int new_carrier_ID;
    //make sure the carrier_reservoir is not empty
    if(state->ReservoirEmpty()){
        state->Grow(eventinfo->growsize, eventinfo->maxpairdegree);
        Grow_non_injection_eventvector(state, eventinfo);
        non_injection_rates->resize(_non_injection_events.size());
    }

    //"buy" the "new" carrier
    new_carrier_ID = state->Buy();
    CarrierBulk* new_carrier = state->GetCarrier(new_carrier_ID);

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
    if(eventinfo->coulomb_strength>0.0) Add_to_mesh(new_carrier_ID,node->position(), eventinfo);

    //Determine effect on Coulomb potential and rates
    Effect_potential_and_rates((int) Add, new_carrier, node, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);
}

void Events:: Remove_carrier(Node* node, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates,
                              Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo) {

    CarrierBulk* removed_carrier = state->GetCarrier(node->occ());

    //remove from graph
    node->RemoveCarrier();    

    //Determine effect on Coulomb potential and rates
    Effect_potential_and_rates((int) Remove, removed_carrier, node, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);    

    
    //remove from mesh
    if(eventinfo->coulomb_strength>0.0) Remove_from_mesh(removed_carrier->id(), node->position(), eventinfo);

    //push back to reservoir
    state->Sell(removed_carrier->id());    
}

inline void Events::Effect_potential_and_rates(int action, CarrierBulk* carrier, Node* node, GraphKMC* graph, StateReservoir* state, Longrange* longrange, 
                                   Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo) {
    Effect_potential_and_non_injection_rates(action, carrier, node, state, longrange, non_injection_rates, eventinfo);

    if(eventinfo->device){
        votca::tools::vec nodeposition = node->position();
        if(nodeposition.z()<eventinfo->coulomb_cut_off_radius) {
            double dist_to_electrode = nodeposition.z();
            Effect_injection_rates(action, carrier, node, graph->left(), dist_to_electrode, state, longrange, left_injection_rates, eventinfo);
        }
        else if(eventinfo->simboxsize.z()-nodeposition.z() < eventinfo->coulomb_cut_off_radius) {
            double dist_to_electrode = eventinfo->simboxsize.z()-nodeposition.z();
            Effect_injection_rates(action, carrier, node, graph->right(), dist_to_electrode, state, longrange, right_injection_rates, eventinfo);
        }
    }

}

void Events::Effect_potential_and_non_injection_rates(int action, CarrierBulk* carrier1, Node* node, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Eventinfo* eventinfo) {

    if(eventinfo->coulomb_strength == 1000.0)
    {
/*        // update event rates for carrier 1 in the case that no coulomb interactions are taken into account
        for (int it = 0; it < node->links().size(); it++) {
            int event_ID1 = carrier1->id()*eventinfo->maxpairdegree+it;
            Node* node2 = node->links()[it]->node2();
            int occ_node = node2->occ();
            bool neighbour = false;
            int event_ID2;
            CarrierBulk* carrier2;
            int reverse_ID;
            if(occ_node != -1 && !eventinfo->no_blocking) {
                carrier2 = state->GetCarrier(occ_node);
                reverse_ID = node->links()[it]->reverse_id();
                event_ID2 = carrier2->id()*eventinfo->maxpairdegree+reverse_ID;
                neighbour = true;
            }
            if(action == (int) Add) {
                _non_injection_events[event_ID1]->Set_event(node->links()[it], carrier1->type(), state, longrange, eventinfo);
                non_injection_rates->setrate(event_ID1, _non_injection_events[event_ID1]->rate());
            }
            else {
                _non_injection_events[event_ID1]->Set_not_in_box_event();
                non_injection_rates->setrate(event_ID1, 0.0);
            }
            if(neighbour && !eventinfo->no_blocking) {
                _non_injection_events[event_ID2]->Set_event(node2->links()[reverse_ID], carrier2->type(), state, longrange, eventinfo);
                non_injection_rates->setrate(event_ID2, _non_injection_events[event_ID2]->rate());                
            }
        }*/
    }
    else {
    
        int interact_sign;
        if(action == (int) Add)    {interact_sign =  1;}
        if(action == (int) Remove) {interact_sign = -1;}

        if(carrier1->type() == (int) Electron) {interact_sign *= -1;}
        if(carrier1->type() == (int) Hole)     {interact_sign *=  1;}

        double RCSQR = eventinfo->coulomb_cut_off_radius*eventinfo->coulomb_cut_off_radius;

        if(eventinfo->device) longrange->Add_charge(1.0*interact_sign, dynamic_cast<NodeDevice*>(node)->layer());

        votca::tools::vec carrier1_pos = node->position();

        // Define cubic boundaries in non-periodic coordinates
        double ix1 = carrier1_pos.x()-eventinfo->coulomb_cut_off_radius-eventinfo->hopdist; double ix2 = carrier1_pos.x()+eventinfo->coulomb_cut_off_radius+eventinfo->hopdist;
        double iy1 = carrier1_pos.y()-eventinfo->coulomb_cut_off_radius-eventinfo->hopdist; double iy2 = carrier1_pos.y()+eventinfo->coulomb_cut_off_radius+eventinfo->hopdist;
        double iz1 = carrier1_pos.z()-eventinfo->coulomb_cut_off_radius-eventinfo->hopdist; double iz2 = carrier1_pos.z()+eventinfo->coulomb_cut_off_radius+eventinfo->hopdist;

        // Break periodicity in x-direction
        if(eventinfo->device) {
            if (iz1<0.0) iz1 = 0.0;
            if (iz2>eventinfo->simboxsize.z()) iz2 = eventinfo->simboxsize.z()-0.5*eventinfo->mindist;
        }
        // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
        int sx1 = floor(ix1/_meshsize_x);
        int sx2 = floor(ix2/_meshsize_x);
        if (ix2 == eventinfo->simboxsize.x()) {sx2 = eventinfo->mesh_size_x;}

        int sy1 = floor(iy1/_meshsize_y);
        int sy2 = floor(iy2/_meshsize_y);
        if (iy2 == eventinfo->simboxsize.y()) {sy2 = eventinfo->mesh_size_y;}

        int sz1 = floor(iz1/_meshsize_z);
        int sz2 = floor(iz2/_meshsize_z);
        if (iz2 == eventinfo->simboxsize.z()) {sz2 = eventinfo->mesh_size_z;}

        // Now visit all relevant sublattices
        for (int isz=sz1; isz<=sz2; isz++) {
            int r_isz = isz;
            while (r_isz < 0) r_isz += eventinfo->mesh_size_z;
            while (r_isz >= eventinfo->mesh_size_z) r_isz -= eventinfo->mesh_size_z;
            for (int isy=sy1; isy<=sy2; isy++) {
                int r_isy = isy;
                while (r_isy < 0) r_isy += eventinfo->mesh_size_y;
                while (r_isy >= eventinfo->mesh_size_y) r_isy -= eventinfo->mesh_size_y;
                for (int isx=sx1; isx<=sx2; isx++) {
                    int r_isx = isx;
                    while (r_isx < 0) r_isx += eventinfo->mesh_size_x;
                    while (r_isx >= eventinfo->mesh_size_x) r_isx -= eventinfo->mesh_size_x;

                    // Ask a list of all charges in this sublattice
                    list<int>::iterator li1,li2,li3;
                    list<int> *carrierlist = &_non_injection_events_mesh[r_isx][r_isy][r_isz];
                    li1 = carrierlist->begin();
                    li2 = carrierlist->end();
                    for (li3=li1; li3!=li2; li3++) {

                        int carrier2_ID = *li3;
                        CarrierBulk* carrier2 = state->GetCarrier(carrier2_ID);

                        int carrier2_type = carrier2->type();
                        Node* carrier2_node = carrier2->node();
                        votca::tools::vec carrier2_pos = carrier2_node->position();

                        int carrier2_charge=0;
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

                        if (carrier2_ID!=carrier1->id()) { // self_image_potential interaction taken into account elsewhere (graph.h)
                            if((node->id()!=carrier2_node->id())&&(distancesqr<=RCSQR)) { 

                                // Charge interacting with its own images, taken care off in graph.h

                                //First we take the direct sr interactions into account
                                if (action == (int) Add) carrier1->Add_from_Coulomb(interact_sign*Compute_Coulomb_potential(np_carrier2_pos.z(),-1.0*distance, true, eventinfo->simboxsize,eventinfo));
                                carrier2->Add_from_Coulomb(interact_sign*Compute_Coulomb_potential(carrier1_pos.z(),distance, true, eventinfo->simboxsize,eventinfo));

                                // This makes sure that charges in the band Rc < x < Rc + hopdist are taken into account correctly (not here)
                            }
                            if (action== (int) Add) {
                                // Adjust Coulomb potential for neighbours of the added carrier
                                for (unsigned it = 0 ; it < node->links().size(); it++) {
                                    votca::tools::vec jumpdistancevector = node->links()[it]->r12();
                                    votca::tools::vec jump_from_carrier1_pos = carrier1_pos + jumpdistancevector;
                                    votca::tools::vec jumpdistance = np_carrier2_pos-jump_from_carrier1_pos;
                                    double distancejumpsqr = abs(jumpdistance)*abs(jumpdistance);

                                    if(distancejumpsqr <= RCSQR) { // makes sure the charge and hopped to charge are in the sphere
                                        if (node->links()[it]->node2()->id() != carrier2_node->id()) { // coulomb interaction on equal nodes lead to zero distance vector (fixed by exciton binding energy)
                                            carrier1->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(np_carrier2_pos.z(),-1.0*jumpdistance, true, eventinfo->simboxsize, eventinfo),it);
                                            //carrier1->Add_to_Coulomb(0.0,it);
                                        }
                                        else {
                                            carrier1->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(np_carrier2_pos.z(),-1.0*jumpdistance, false, eventinfo->simboxsize, eventinfo),it);
                                            //carrier1->Add_to_Coulomb(0.0,it);
                                        }
                                    }
                                }
                            }
                            // Adjust Coulomb potential and event rates for neighbours of carrier2
                            std::vector<Link*>::iterator it;
                            for (unsigned it = 0; it < carrier2_node->links().size(); it++) {
                                votca::tools::vec jumpdistancevector = carrier2_node->links()[it]->r12();
                                votca::tools::vec jump_from_carrier2_pos = np_carrier2_pos+jumpdistancevector;
                                votca::tools::vec jumpdistance = jump_from_carrier2_pos - carrier1_pos;
                                double distancejumpsqr = abs(jumpdistance)*abs(jumpdistance);
                                int event_ID = carrier2->id()*eventinfo->maxpairdegree + it;

                                if(distancejumpsqr <= RCSQR) {
                                    if(carrier2_node->links()[it]->node2()->id() == node->id()) {
                                        carrier2->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(carrier1_pos.z(),jumpdistance,false,eventinfo->simboxsize, eventinfo), it);
                                        //carrier2->Add_to_Coulomb(0.0, it);
                                        _non_injection_events[event_ID]->Set_event(carrier2_node->links()[it], carrier2_type, state, longrange, eventinfo); // if carrier2 is actually linking with carrier1
                                    }
                                    else {
                                        carrier2->Add_to_Coulomb(interact_sign*Compute_Coulomb_potential(carrier1_pos.z(),jumpdistance,true,eventinfo->simboxsize, eventinfo), it);
                                        //carrier2->Add_to_Coulomb(0.0, it);
                                        _non_injection_events[event_ID]->Determine_rate(state, longrange, eventinfo);
                                    }
                                }
                                non_injection_rates->setrate(event_ID,_non_injection_events[event_ID]->rate());
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

          // update event rates for carrier 1 , done after all carriers within radius coulomb_cut_off_radius are checked
        for (unsigned it = 0; it < node->links().size(); it++) {
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

}

void Events:: Effect_injection_rates(int action, CarrierBulk* carrier, Node* node, Node* electrode, double dist_to_electrode, StateReservoir* state, Longrange* longrange, Bsumtree* injection_rates, Eventinfo* eventinfo){

    int interact_sign;
    //int z_mesh;
    
    if(action == (int) Add) {interact_sign = 1;}
    if(action == (int) Remove) {interact_sign = -1;}
    if(carrier->type() == (int) Electron) {interact_sign *= -1;}
    if(carrier->type() == (int) Hole) {interact_sign *= 1;}
    
    votca::tools::vec carrier1_pos = node->position();
//    votca::tools::vec flipped_pos;
//    if(electrode->type() == LeftElectrodeNode) {flipped_pos = carrier1_pos;}
//    else {flipped_pos = votca::tools::vec(eventinfo->simboxsize.x()-carrier1_pos.x(),carrier1_pos.y(),carrier1_pos.z());}
    double RCSQR = eventinfo->coulomb_cut_off_radius*eventinfo->coulomb_cut_off_radius;

    // typically one wants hopdist to be smaller than RC
    double bound = eventinfo->coulomb_cut_off_radius;
//    if (flipped_pos.x()>eventinfo->hopdist)       bound = sqrt(double(RCSQR - (dist_to_electrode-eventinfo->hopdist)*(dist_to_electrode-eventinfo->hopdist)));
//    else if (flipped_pos.x()<=eventinfo->hopdist) bound = eventinfo->coulomb_cut_off_radius;

    // Define cubic boundaries in non-periodic coordinates
    double ix1 = carrier1_pos.x()-bound-eventinfo->hopdist; double ix2 = carrier1_pos.x()+bound+eventinfo->hopdist;
    double iy1 = carrier1_pos.y()-bound-eventinfo->hopdist; double iy2 = carrier1_pos.y()+bound+eventinfo->hopdist;
    double iz1 = 0.0;                                       double iz2 = eventinfo->hopdist;
    
    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    
    int sx1 = floor(ix1/_meshsize_x); int sx2 = floor(ix2/_meshsize_x);
    int sy1 = floor(iy1/_meshsize_y); int sy2 = floor(iy2/_meshsize_y);
    int sz1 = floor(iz1/_meshsize_z); int sz2 = floor(iz2/_meshsize_z);

    if (ix2 == eventinfo->simboxsize.x()) {sx2 = eventinfo->mesh_size_x;}    
    if (iy2 == eventinfo->simboxsize.y()) {sy2 = eventinfo->mesh_size_y;}    

    
//    if(sy2 == eventinfo->mesh_size_y) {sy2--;}
//    if(sz2 == eventinfo->mesh_size_z) {sz2--;}    
   
    // Now visit all relevant sublattices
    for (int isx=sx1; isx<=sx2; isx++) {
        int r_isx = isx;
        while (r_isx < 0) r_isx += eventinfo->mesh_size_x;
        while (r_isx >= eventinfo->mesh_size_x) r_isx -= eventinfo->mesh_size_x;
        for (int isy=sy1; isy<=sy2; isy++) {
            int r_isy = isy;
            while (r_isy < 0) r_isy += eventinfo->mesh_size_y;
            while (r_isy >= eventinfo->mesh_size_y) r_isy -= eventinfo->mesh_size_y;
            for (int isz = sz1; isz<=sz2; isz++) {
                int r_isz = isz;

                // Ask a list of all nodes in this sublattice
                list<int>::iterator li1,li2,li3;
                list<int> *eventmesh ;
                
                int IDelectrodeswitch;

                if(electrode->type() == (int) LeftElectrodeNode) {
                    eventmesh = &_left_injection_events_mesh[r_isx][r_isy][r_isz];
                    IDelectrodeswitch = 0;
                }
                else if(electrode->type() == (int) RightElectrodeNode) {
                    eventmesh = &_right_injection_events_mesh[r_isx][r_isy][r_isz];
                    IDelectrodeswitch = _total_left_injection_events;
                }
                else{
                    eventmesh=NULL;
                }

                li1 = eventmesh->begin();
                li2 = eventmesh->end();
                for (li3=li1; li3!=li2; li3++) {
                    int eventID = *li3;
                    Event* inject_event = _injection_events[eventID];
                    Node* probenode = inject_event->link()->node2();
                    votca::tools::vec probepos = probenode->position();
         
                    // Compute coordinates in non-periodic lattice
          
                    votca::tools::vec periodic_convert = votca::tools::vec((isx-r_isx)*_meshsize_x,(isy-r_isy)*_meshsize_y,0.0);
                    votca::tools::vec np_probepos = probepos + periodic_convert;
                    votca::tools::vec distance = np_probepos-carrier1_pos;

                    double distancesqr = abs(distance)*abs(distance);
                    if(distancesqr <= RCSQR) {

                        if(inject_event->carrier_type() == (int) Electron) {interact_sign *= -1;}
                    
                        if (probenode->id()!=node->id()) {
                            inject_event->Add_injection_potential(interact_sign*Compute_Coulomb_potential(carrier1_pos.z(),distance,true,eventinfo->simboxsize,eventinfo));                    
                            //inject_event->Add_injection_potential(0.0);                    
                            _injection_events[eventID]->Determine_rate(state, longrange, eventinfo);
                            injection_rates->setrate(eventID - IDelectrodeswitch, _injection_events[eventID]->rate());
                        }
                        else {
                            inject_event->Add_injection_potential(interact_sign*Compute_Coulomb_potential(carrier1_pos.z(),distance,false,eventinfo->simboxsize,eventinfo));                    
                            //inject_event->Add_injection_potential(0.0);                    
                            _injection_events[eventID]->Set_event(inject_event->link(), inject_event->carrier_type(),state, longrange, eventinfo);
                            injection_rates->setrate(eventID - IDelectrodeswitch, _injection_events[eventID]->rate());
                        }
                    }
                }
            }
        }
    }  
}

double Events::Compute_Coulomb_potential(double startz, votca::tools::vec dif, bool direct, votca::tools::vec simboxsize, Eventinfo* eventinfo) {

    double RC = eventinfo->coulomb_cut_off_radius;
    double RCSQR = RC*RC;
    double coulpot;
    
    if(direct){
        coulpot = 1.0/abs(dif) - 1.0/sqrt(RCSQR);
    }
    else {
        coulpot = 0.0;
    }
    
    if(eventinfo->device) {
        
        double L = simboxsize.z();
        double distsqr_planar = dif.x()*dif.x() + dif.y()*dif.y();
      
        double sign;
        double distz_1;
        double distz_2;
        double distancesqr_1;
        double distancesqr_2;
       // bool outside_cut_off1 = false;
        //bool outside_cut_off2 = false;
      
        for (int i=0;i<eventinfo->number_short_range_images; i++) {
            if (div(i,2).rem==0) { // even generation
                sign = -1.0;
                distz_1 = i*L + 2*startz + dif.z();
                distz_2 = (i+2)*L - 2*startz - dif.z(); 
            }
            else {
                sign = 1.0;
                distz_1 = (i+1)*L + dif.z();
                distz_2 = (i+1)*L - dif.z();
            }
            distancesqr_1 = distz_1*distz_1 + distsqr_planar;
            if(distancesqr_1 <= RCSQR) coulpot += sign*(1.0/sqrt(distancesqr_1) - 1.0/sqrt(RCSQR));
            distancesqr_2 = distz_2*distz_2 + distsqr_planar;
            if(distancesqr_2 <= RCSQR) coulpot += sign*(1.0/sqrt(distancesqr_2) - 1.0/sqrt(RCSQR));
        }

    }
    return coulpot;
}

void Events::Recompute_all_events(StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo) {
  
    Recompute_all_non_injection_events(state,longrange, non_injection_rates, eventinfo);
    Recompute_all_injection_events(state,longrange,left_injection_rates, right_injection_rates, eventinfo);
}

void Events::Recompute_all_non_injection_events(StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Eventinfo* eventinfo) {
    
    std::vector<Event*>::iterator it;
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

void Events::Recompute_all_injection_events(StateReservoir* state, Longrange* longrange, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo) {
    
    std::vector<Event*>::iterator it;
    for (it = _injection_events.begin(); it!=_injection_events.end(); it++){
        (*it)->Determine_rate(state, longrange, eventinfo);
        if((*it)->link()->node1()->type() == (int) LeftElectrodeNode) {
            left_injection_rates->setrate((*it)->id(),(*it)->rate());       
        } 
        else {
            right_injection_rates->setrate((*it)->id() - _total_left_injection_events,(*it)->rate()); 
        }
    }
}


void Events::Initialize_device_eventvector(GraphKMC* graph, StateReservoir* state, Longrange* longrange, Eventinfo* eventinfo){ //

    this->Initialize_bulk_eventvector(graph, state, eventinfo);
    
    _injection_events.clear();
    int Event_id_count = 0;
    if(eventinfo->left_electron_injection) {
        Initialize_injection_eventvector(Event_id_count,graph->left(), (int) Electron, state, longrange, eventinfo); 
        Event_id_count += graph->left()->links().size();
        _total_left_injection_events += graph->left()->links().size();
    }
    if(eventinfo->left_hole_injection) {
        Initialize_injection_eventvector(Event_id_count,graph->left(), (int) Hole, state, longrange, eventinfo); 
        Event_id_count += graph->left()->links().size();
        _total_left_injection_events += graph->left()->links().size();
    }
    if(eventinfo->right_electron_injection) {
        Initialize_injection_eventvector(Event_id_count,graph->right(), (int) Electron, state, longrange, eventinfo); 
        Event_id_count += graph->right()->links().size();
        _total_right_injection_events += graph->right()->links().size();

    }
    if(eventinfo->right_hole_injection) {
        Initialize_injection_eventvector(Event_id_count,graph->right(), (int) Hole, state, longrange, eventinfo);
        _total_right_injection_events += graph->right()->links().size();
    }
}

void Events::Initialize_bulk_eventvector(GraphKMC* graph, StateReservoir* state, Eventinfo* eventinfo)
{
    std::vector<Event*>::iterator it;

    _non_injection_events.clear();
    Grow_non_injection_eventvector(state, eventinfo);

    _total_left_injection_events = 0;
    _total_right_injection_events = 0;    
    _total_non_injection_events = 0;    
}

void Events::Initialize_rates(Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo){

    non_injection_rates->initialize(_non_injection_events.size());
    std::vector<Event*>::iterator it; 
    for (it = _non_injection_events.begin(); it!=_non_injection_events.end(); it++) {non_injection_rates->setrate((*it)->id(),(*it)->rate());}
  
    if(eventinfo->device) {
        left_injection_rates->initialize(_total_left_injection_events);
        right_injection_rates->initialize(_total_right_injection_events);

        for (it = _injection_events.begin(); it!=_injection_events.begin()+_total_left_injection_events; it++) {left_injection_rates->setrate((*it)->id(),(*it)->rate());}  
        for (it = _injection_events.begin()+_total_left_injection_events; it!=_injection_events.end(); it++) {right_injection_rates->setrate((*it)->id() - _total_left_injection_events,(*it)->rate());}
    }    
}

void Events::Initialize_injection_eventvector(int Event_id_count, Node* electrode, int carrier_type, StateReservoir* state, Longrange* longrange, Eventinfo* eventinfo){

    int Event_map = Event_id_count;
    
    for (unsigned it = 0; it < electrode->links().size(); it++) {

        Event *newEvent = new Event(Event_map, electrode->links()[it], carrier_type, state, longrange, eventinfo);
        _injection_events.push_back(newEvent);
        Event_map++;
        
    }
}

void Events::Initialize_after_charge_placement(GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo){

    for(int icar = 0; icar< state->GetCarrierSize(); icar++) {
        CarrierBulk* probe_carrier = state->GetCarrier(icar);
        Node* carrier_node = probe_carrier->node();
        if(probe_carrier->inbox()) {
           if(eventinfo->coulomb_strength>0.0) Add_to_mesh(icar,carrier_node->position(), eventinfo);
           Effect_potential_and_rates((int) Add, probe_carrier, carrier_node, graph, state, longrange, non_injection_rates, left_injection_rates, right_injection_rates, eventinfo);
        }
    }
}


void Events::Grow_non_injection_eventvector(StateReservoir* state, Eventinfo* eventinfo){

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
        //CarrierBulk* probecarrier = state->GetCarrier(carrier_ID);
        int Event_map = carrier_ID*eventinfo->maxpairdegree;
        Event *newEvent;

        for(int ifill = 0; ifill<eventinfo->maxpairdegree; ifill++) {
            newEvent = new Event(Event_map, (int) Notinbox); 
            _non_injection_events.push_back(newEvent); 
            Event_map ++; 
        } // carrier not in box
    }
}

void Events::Init_non_injection_meshes(Eventinfo* eventinfo) {

    // determine the dimensions of the meshes
  
    votca::tools::vec simboxsize = eventinfo->simboxsize;
    _meshsize_x = simboxsize.x()/eventinfo->mesh_size_x; 
    _meshsize_y = simboxsize.y()/eventinfo->mesh_size_y; 
    _meshsize_z = simboxsize.z()/eventinfo->mesh_size_z;

    _non_injection_events_mesh = Resize_mesh(eventinfo->mesh_size_x,eventinfo->mesh_size_y,eventinfo->mesh_size_z);

}

void Events::Init_injection_meshes(StateReservoir* state, Eventinfo* eventinfo) {

    _inject_meshnr_z = floor(eventinfo->hopdist/_meshsize_z)+1;

    _left_injection_events_mesh = Resize_mesh(eventinfo->mesh_size_x,eventinfo->mesh_size_y,_inject_meshnr_z);
    _right_injection_events_mesh = Resize_mesh(eventinfo->mesh_size_x,eventinfo->mesh_size_y,_inject_meshnr_z);    
    std::vector<Event*>::iterator it;
    for (it = _injection_events.begin(); it != _injection_events.end(); it++ ) {
        if((*it)->link()->node1()->type()==LeftElectrodeNode) {
            votca::tools::vec position = (*it)->link()->node2()->position();
            
            double posx = position.x(); int iposx = floor(posx/_meshsize_x);
            double posy = position.y(); int iposy = floor(posy/_meshsize_y);
            double posz = position.z(); int iposz = floor(posz/_meshsize_z);

            if(iposx == eventinfo->mesh_size_x) {iposx--;}
            if(iposy == eventinfo->mesh_size_y) {iposy--;}           
            _left_injection_events_mesh[iposx][iposy][iposz].push_back((*it)->id());
        }
        
        //for the righthandside electrode, we take the distance from said electrode

        if((*it)->link()->node1()->type()==RightElectrodeNode) {

            votca::tools::vec simboxsize = eventinfo->simboxsize;
            votca::tools::vec unflipped_position = (*it)->link()->node2()->position();
            votca::tools::vec position = votca::tools::vec(unflipped_position.x(),unflipped_position.y(),simboxsize.z() -unflipped_position.z());

            double posx = position.x(); int iposx = floor(posx/_meshsize_x);
            double posy = position.y(); int iposy = floor(posy/_meshsize_y);
            double posz = position.z(); int iposz = floor(posz/_meshsize_z);

            if(iposx == eventinfo->mesh_size_x) {iposx--;}
            if(iposy == eventinfo->mesh_size_y) {iposy--;}
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

    if(iposx == eventinfo->mesh_size_x) {iposx--;}
    if(iposy == eventinfo->mesh_size_y) {iposy--;}
    if(iposz == eventinfo->mesh_size_z) {iposz--;}
    
    _non_injection_events_mesh[iposx][iposy][iposz].push_back(ID);
};

void Events::Remove_from_mesh(int ID, votca::tools::vec position, Eventinfo* eventinfo){
    
    double posx = position.x(); int iposx = floor(posx/_meshsize_x);
    double posy = position.y(); int iposy = floor(posy/_meshsize_y);
    double posz = position.z(); int iposz = floor(posz/_meshsize_z);

    if(iposx == eventinfo->mesh_size_x) {iposx--;}
    if(iposy == eventinfo->mesh_size_y) {iposy--;}
    if(iposz == eventinfo->mesh_size_z) {iposz--;}
         
    _non_injection_events_mesh[iposx][iposy][iposz].remove(ID);
}

double Events::av_inject_energyfactor(int electrode) {
    double av_energyfactor = 0.0;
    if(electrode == 0) {
        for(int ievent = 0; ievent < _total_left_injection_events ; ievent++) {
            av_energyfactor += _injection_events[ievent]->energyfactor();
        }
    }
    else {
        for(int ievent = _total_left_injection_events; ievent < _total_left_injection_events + _total_right_injection_events ; ievent++) {
            av_energyfactor += _injection_events[ievent]->energyfactor();
        }        
    }
    
    return av_energyfactor;
}

double Events::av_rate(int electrode) {
    double av_rate = 0.0;
    if(electrode == 0) {
        for(int ievent = 0; ievent < _total_left_injection_events ; ievent++) {
            av_rate += _injection_events[ievent]->rate();
        }
        av_rate /= _total_left_injection_events;
    }
    else {
        for(int ievent = _total_left_injection_events; ievent < _total_left_injection_events + _total_right_injection_events ; ievent++) {
            av_rate += _injection_events[ievent]->rate();
        }        
        av_rate /= _total_right_injection_events;
    }
    
    return av_rate;
}

double Events::av_inject_transferfactor(int electrode) {
    double av_transferfactor = 0.0;
    if(electrode == 0) {
        for(int ievent = 0; ievent < _total_left_injection_events ; ievent++) {
            av_transferfactor += _injection_events[ievent]->transferfactor();
        }
    }
    else {
        for(int ievent = _total_left_injection_events; ievent < _total_left_injection_events + _total_right_injection_events ; ievent++) {
            av_transferfactor += _injection_events[ievent]->transferfactor();
        }        
    }
    return av_transferfactor;
}

        
    }
}
