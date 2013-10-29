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
    State* state;
    Graph* graph;
    Longrange* longrange;

    void Effect_potential_and_non_injection_rates_one_charge(action AR, Carrier* carrier,
                                                   int meshsizeX, int meshsizeY, int meshsizeZ, vector< vector< vector< vector <list<int> > > > > coulomb_mesh,
                                                   vector<Carrier*> electrons, vector<Carrier*> holes, double coulcut, double hopdist, bool device, myvec sim_box_size,
                                                   int nr_sr_images, vector<Node*> nodes, int max_pair_degree, string formalism, Globaleventinfo* globinfo);
    void Effect_injection_rates_one_charge(action AR, vector<Node*> nodes, Carrier* carrier, 
                                                   vector< vector< vector <list<Node*> > > > node_mesh,
                                                   int meshsizeX, int meshsizeY, int meshsizeZ,
                                                   double coulcut, double hopdist, double dist_to_electrode, Node* electrode, 
                                                   myvec sim_box_size, int nr_sr_images, string formalism, Globaleventinfo* globinfo, bool dual_injection, int nr_left_injector_nodes);
    
    
    void Recompute_all_injection_events(bool dual_injection, Node* left_electrode, Node* right_electrode, int nr_left_injector_nodes, string formalism, Globaleventinfo* globevent);
    void Recompute_all_non_injection_events(vector<Node*> nodes, vector<Carrier*> electrons, vector<Carrier*> holes, int maxpairdegree,  string formalism, Globaleventinfo* globevent, bool device);
  
    void Initialize_eventvector_for_device(bool dual_injection, vector <Carrier*> electrons, vector <Carrier*> holes, Node* left_electrode, Node* right_electrode, int maxpairdegree); // if dual_injection is true, both types of carriers are injected from both electrodes
    void Grow_non_injection_eventvector(int carrier_grow_size, vector<Carrier*> carriers, vector<Event*> eventvector,int maxpairdegree);
    
    double Compute_Coulomb_potential(double startx, myvec dif, myvec sim_box_size, double coulcut, int nr_sr_images, bool device);
    
private:
    void Initialize_injection_eventvector(Node* electrode, vector<Event*> eventvector);
    
};


void Events::Effect_potential_and_non_injection_rates_one_charge(action AR, Carrier* carrier,
                                                   int meshsizeX, int meshsizeY, int meshsizeZ, vector< vector< vector< vector <list<int> > > > > coulomb_mesh,
                                                   vector<Carrier*> electrons, vector<Carrier*> holes, double coulcut, double hopdist, bool device, myvec sim_box_size,
                                                   int nr_sr_images, vector<Node*> nodes, int max_pair_degree, string formalism, Globaleventinfo* globevent) {

    int interact_sign;
    
    if(AR == Add) {interact_sign = 1;}
    if(AR == Remove) {interact_sign = -1;}
    if(carrier->carrier_type == Electron) {interact_sign *= -1;}
    if(carrier->carrier_type == Hole) {interact_sign *=1;}
    
    Node* carnode = graph->nodes[carrier->carrier_node_ID];
    
    //calculate the change to the longrange cache
    
    if(device){ 
        int layer_index = carnode->layer_index;
        longrange->layercharge[layer_index] += interact_sign;
    }
     
    myvec carpos = carnode->node_position;

    // Define cubic boundaries in non-periodic coordinates
    double ix1 = carpos.x()-coulcut-hopdist; double ix2 = carpos.x()+coulcut+hopdist;
    double iy1 = carpos.y()-coulcut-hopdist; double iy2 = carpos.y()+coulcut+hopdist;
    double iz1 = carpos.z()-coulcut-hopdist; double iz2 = carpos.z()+coulcut+hopdist;

    // Break periodicity in x-direction
    if(device) {
        if (ix1<0.0) ix1 = 0.0;
        if (ix2>=sim_box_size.x()) ix2 = sim_box_size.x();
    }
  
    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    int sx1 = floor(ix1/coulcut);
    int sx2 = floor(ix2/coulcut);
    int sy1 = floor(iy1/coulcut);
    int sy2 = floor(iy2/coulcut);
    int sz1 = floor(iz1/coulcut);
    int sz2 = floor(iz2/coulcut);
  
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
        int r_isz = isz;
        while (r_isz < 0) r_isz += meshsizeZ;
        while (r_isz >= meshsizeZ) r_isz -= meshsizeZ;
        for (int isy=sy1; isy<=sy2; isy++) {
            int r_isy = isy;
            while (r_isy < 0) r_isy += meshsizeY;
            while (r_isy >= meshsizeY) r_isy -= meshsizeY;
            for (int isx=sx1; isx<=sx2; isx++) {
                int r_isx = isx;
                while (r_isx < 0) r_isx += meshsizeX;
                while (r_isx >= meshsizeX) r_isx -= meshsizeX;
                for (int icartype = 0;icartype <2;icartype++) {
                
                    // Ask a list of all charges in this sublattice
                    list<int>::iterator li1,li2,li3;
                    list<int> *carrierList = &coulomb_mesh[r_isx][r_isy][r_isz][icartype];
                    li1 = carrierList->begin();
                    li2 = carrierList->end();
                    for (li3=li1; li3!=li2; li3++) {
                        int probecarrier_ID = *li3;
                        Carrier* probecarrier;
                        if(icartype == 0) {probecarrier = electrons[probecarrier_ID];}
                        if(icartype == 1) {probecarrier = holes[probecarrier_ID];}
                        Node* probenode = nodes[probecarrier->carrier_node_ID];
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
                        myvec periodic_convert = myvec((isx-r_isx)*coulcut,(isy-r_isy)*coulcut,(isz-r_isz)*coulcut);
                        myvec np_probepos = probepos + periodic_convert;
                        myvec distance = np_probepos-carpos;

                        double distancesqr = abs(distance)*abs(distance);

                        if (probecarrier_ID!=carrier->carrier_ID) {
                            if((carnode->node_ID!=probenode->node_ID)&&(distancesqr<=coulcut*coulcut)) { 
                                
                                // Charge interacting with its own images, taken care off in graph.h
                                // In case multiple charges are on the same node, coulomb calculation on the same spot is catched
                          
                                //First we take the direction sr interactions into account
                                if (AR==Add) carrier->srfrom +=interact_sign*Compute_Coulomb_potential(np_probepos.x(),-1.0*distance,
                                                            sim_box_size, coulcut, nr_sr_images, device);
                                probecarrier->srfrom += interact_sign*Compute_Coulomb_potential(carpos.x(),distance,
                                                            sim_box_size, coulcut, nr_sr_images, device);
                            }
                            if (AR==Add) {
              
                                // Adjust Coulomb potential for neighbours of the added carrier
                                for (int jump=0; jump < carnode->pairing_nodes.size(); jump++) {
                                    myvec jumpdistancevector = carnode->static_event_info[jump].distance;
                                    myvec jumpcarrierpos = carnode->node_position + jumpdistancevector;
                                    myvec jumpdistance = np_probepos - jumpcarrierpos;
                                    double distancejumpsqr = abs(jumpdistance)*abs(jumpdistance);

                                    if(distancejumpsqr <= coulcut*coulcut) {
                                    
                                        carrier->srto[jump] += interact_sign*Compute_Coulomb_potential(np_probepos.x(),jumpdistance,
                                                         sim_box_size, coulcut, nr_sr_images, device);
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
           
                            // Adjust Coulomb potential for neighbours of carrier2
                            for (int jump=0; jump < probenode->pairing_nodes.size(); jump++) {
                                myvec jumpdistancevector = probenode->static_event_info[jump].distance;
                                myvec jumpprobepos = np_probepos+jumpdistancevector;
                                myvec jumpdistance = carpos-jumpprobepos;
                                double distsqr = abs(jumpdistance)*abs(jumpdistance);
                                int event_ID = probecarrier->carrier_ID*max_pair_degree+jump;
                                
                                double fromlongrange;
                                double tolongrange;
                                if(device) {
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
                                
                                if(distsqr <= coulcut*coulcut) {
                                    if(probecarrier->carrier_type==Electron) {
                                        probecarrier->srto[jump] += 
                                                        interact_sign*Compute_Coulomb_potential(carpos.x(),jumpdistance,
                                                        sim_box_size, coulcut, nr_sr_images, device);
                                        
                                        El_non_injection_events[event_ID]->Set_non_injection_event(nodes, probecarrier, jump, formalism, fromlongrange, tolongrange, globevent);
                                        El_non_injection_rates->setrate(event_ID, El_non_injection_events[event_ID]->rate);
                                    }
                                    else if(probecarrier->carrier_type==Hole) {
                                        probecarrier->srto[jump] += 
                                                            interact_sign*Compute_Coulomb_potential(carpos.x(),jumpdistance,
                                                            sim_box_size, coulcut, nr_sr_images, device);
                                        Ho_non_injection_events[event_ID]->Set_non_injection_event(nodes, probecarrier, jump, formalism, fromlongrange, tolongrange, globevent);
                                        Ho_non_injection_rates->setrate(event_ID, Ho_non_injection_events[event_ID]->rate);
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
    for (int jump=0; jump < carnode->pairing_nodes.size(); jump++) {
        int event_ID = carrier->carrier_ID*max_pair_degree+jump;
    
        double fromlongrange;
        double tolongrange;
        if(device) {
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
                El_non_injection_events[event_ID]->Set_non_injection_event(nodes, carrier, jump, formalism, fromlongrange, tolongrange, globevent);
                El_non_injection_rates->setrate(event_ID, El_non_injection_events[event_ID]->rate);
            }
            else {
                El_non_injection_events[event_ID]->fromtype = Fromnotinbox;
                El_non_injection_events[event_ID]->totype = Tonotinbox;
                El_non_injection_events[event_ID]->rate = 0.0;
                El_non_injection_rates->setrate(event_ID, 0.0);
            }
        }
        else if(carrier->carrier_type==Hole) {
            if(AR == Add) {
                Ho_non_injection_events[event_ID]->Set_non_injection_event(nodes, carrier, jump, formalism, fromlongrange, tolongrange, globevent);
                Ho_non_injection_rates->setrate(event_ID, Ho_non_injection_events[event_ID]->rate);
            }
            else {
                Ho_non_injection_events[event_ID]->fromtype = Fromnotinbox;
                Ho_non_injection_events[event_ID]->totype = Tonotinbox;
                Ho_non_injection_events[event_ID]->rate = 0.0;
                Ho_non_injection_rates->setrate(event_ID, 0.0);
            }            
        }
    }
}        
        
void Events::Effect_injection_rates_one_charge(action AR, vector<Node*> nodes, Carrier* carrier, 
                                                   vector< vector< vector <list<Node*> > > > node_mesh,
                                                   int meshsizeX, int meshsizeY, int meshsizeZ,
                                                   double coulcut, double hopdist, double dist_to_electrode, Node* electrode, 
                                                   myvec sim_box_size, int nr_sr_images, string formalism, Globaleventinfo* globevent, bool dual_injection, int nr_left_injector_nodes) {
                                                   
    int interact_sign;
    int x_mesh;
    
    if(AR == Add) {interact_sign = 1;}
    if(AR == Remove) {interact_sign = -1;}
    if(carrier->carrier_type == Electron) {interact_sign *= -1;}
    if(carrier->carrier_type == Hole) {interact_sign *= 1;}
    if(electrode->node_type == LeftElectrode) {x_mesh = 0;}
    if(electrode->node_type == RightElectrode) {x_mesh = meshsizeX-1;}
    
    Node* carnode = nodes[carrier->carrier_node_ID];
    myvec carpos = carnode->node_position;
  
    if (dist_to_electrode <= coulcut) { // Distance to electrode
  
        //maximum entry inside device is hopdist (make this variable!!)

        double bound = sqrt(double(coulcut*coulcut - dist_to_electrode*dist_to_electrode));

        // Define cubic boundaries in non-periodic coordinates
        double iy1 = carpos.y()-bound-hopdist; double iy2 = carpos.y()+bound+hopdist;
        double iz1 = carpos.z()-bound-hopdist; double iz2 = carpos.z()+bound+hopdist;

        // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
        int sy1 = floor(iy1/hopdist);
        int sy2 = floor(iy2/hopdist);
        int sz1 = floor(iz1/hopdist);
        int sz2 = floor(iz2/hopdist);
  
        // Now visit all relevant sublattices
        for (int isz=sz1; isz<=sz2; isz++) {
            int r_isz = isz;
            while (r_isz < 0) r_isz += meshsizeZ;
            while (r_isz >= meshsizeZ) r_isz -= meshsizeZ;
            for (int isy=sy1; isy<=sy2; isy++) {
                int r_isy = isy;
                while (r_isy < 0) r_isy += meshsizeY;
                while (r_isy >= meshsizeY) r_isy -= meshsizeY;
       
                // Ask a list of all nodes in this sublattice
                list<Node*>::iterator li1,li2,li3;
                list<Node*> *nodemesh = &node_mesh[x_mesh][r_isy][r_isz];
                li1 = nodemesh->begin();
                li2 = nodemesh->end();
                for (li3=li1; li3!=li2; li3++) {
                    Node* probenode = *li3;
                    myvec probepos = probenode->node_position;
          
                    // Compute coordinates in non-periodic lattice
          
                    myvec periodic_convert = myvec(0.0,(isy-r_isy)*hopdist,(isz-r_isz)*hopdist);
                    myvec np_probepos = probepos + periodic_convert;
                    myvec distance = np_probepos-carpos;

                    double distancesqr = abs(distance)*abs(distance);

                    if ((probenode->node_ID!=carnode->node_ID)&&(distancesqr <= coulcut*coulcut)) { // calculated for holes, multiply interact_sign with -1 for electrons
                        probenode->injection_potential +=interact_sign*Compute_Coulomb_potential(carpos.x(),distance,sim_box_size,coulcut,nr_sr_images,true);
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
                            Ho_injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Hole, formalism, 0.0, tolongrange, globevent);
                            Ho_injection_rates->setrate(event_ID, Ho_injection_events[event_ID]->rate);
                            if(dual_injection) {
                                El_injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Electron, formalism, 0.0, tolongrange, globevent);
                                El_injection_rates->setrate(event_ID, El_injection_events[event_ID]->rate);
                            }
                        }
                        else if(electrode->node_type == RightElectrode) {
                            injector_ID = probenode->right_injector_ID;
                            event_ID = injector_ID;
                            if(dual_injection) {
                                event_ID += nr_left_injector_nodes;
                                Ho_injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Hole, formalism, 0.0, tolongrange, globevent);
                                Ho_injection_rates->setrate(event_ID, Ho_injection_events[event_ID]->rate); 
                            }                               
                            El_injection_events[event_ID]->Set_injection_event(electrode, injector_ID, 
                                                  Electron, formalism, 0.0, tolongrange, globevent);
                            El_injection_rates->setrate(event_ID, El_injection_events[event_ID]->rate);
                        }
                    }
                }
            }
        }
    }  
}

double Events::Compute_Coulomb_potential(double startx, myvec dif, myvec sim_box_size, double coulcut, int nr_sr_images, bool device) {

    double coulpot;
    double RC = coulcut;
    double RCSQR = RC*RC;
   
    if(!device) {
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
            for (int i=0;i<nr_sr_images; i++) {
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


void Events::Recompute_all_non_injection_events(vector<Node*> nodes, vector<Carrier*> electrons, vector<Carrier*> holes, int maxpairdegree, string formalism, Globaleventinfo* globevent, bool device) {
    
    int Event_map;
       
    for (int electron_ID = 0; electron_ID<electrons.size(); electron_ID++) {
        Node* electron_node = nodes[electrons[electron_ID]->carrier_node_ID];
        for (int ipair = 0; ipair < electron_node->pairing_nodes.size();ipair++){
            
            Event_map = electron_ID*maxpairdegree + ipair;
            Carrier* electron = electrons[electron_ID];
            
            double lrfrom;
            double lrto;
            
            if(device && electron->is_in_sim_box){
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
            
            El_non_injection_events[Event_map]->Set_non_injection_event(nodes,electron, ipair, formalism,lrfrom,lrto, globevent);
            El_non_injection_rates->setrate(Event_map,El_non_injection_events[Event_map]->rate);
        }
    }

    for (int hole_ID = 0; hole_ID<holes.size(); hole_ID++) {
        Node* hole_node = nodes[holes[hole_ID]->carrier_node_ID];
        for (int ipair = 0; ipair < hole_node->pairing_nodes.size();ipair++){
            
            Event_map = hole_ID*maxpairdegree + ipair;
            Carrier* hole = holes[hole_ID];
            
            double lrfrom;
            double lrto;
            
            if(device && hole->is_in_sim_box){
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
            
            Ho_non_injection_events[Event_map]->Set_non_injection_event(nodes,hole, ipair, formalism,lrfrom ,lrto, globevent);
            Ho_non_injection_rates->setrate(Event_map,Ho_non_injection_events[Event_map]->rate);
        }
    }
}

void Events::Recompute_all_injection_events(bool dual_injection, Node* left_electrode, Node* right_electrode, int nr_left_injector_nodes, string formalism, Globaleventinfo* globevent) {
    
    int Event_map;

    for (int inject_node = 0; inject_node<left_electrode->pairing_nodes.size(); inject_node++) {

        Event_map = inject_node;

        double lrto;
        if(left_electrode->pairing_nodes[inject_node]->node_type == Normal) {
            lrto = longrange->Get_cached_longrange(left_electrode->pairing_nodes[inject_node]->layer_index);
        }
        else { // Collection
            lrto = 0.0;
        }
        
        El_injection_events[Event_map]->Set_injection_event(left_electrode, inject_node, Electron, formalism, 0.0, lrto, globevent);   
        El_injection_rates->setrate(Event_map,El_injection_events[Event_map]->rate);
        if(dual_injection) {
            Ho_injection_events[Event_map]->Set_injection_event(left_electrode, inject_node, Hole, formalism, 0.0, lrto, globevent);   
            Ho_injection_rates->setrate(Event_map,Ho_injection_events[Event_map]->rate);
        }        
    }
    
    for (int inject_node = 0; inject_node<right_electrode->pairing_nodes.size(); inject_node++) {

        double lrto;
        if(right_electrode->pairing_nodes[inject_node]->node_type == Normal) {
            lrto = longrange->Get_cached_longrange(right_electrode->pairing_nodes[inject_node]->layer_index);
        }
        else { // Collection
            lrto = 0.0;
        }        
        
        
        if(dual_injection){
            Event_map = nr_left_injector_nodes + inject_node;
            El_injection_events[Event_map]->Set_injection_event(right_electrode, inject_node, Electron, formalism, 0.0 , lrto, globevent);
            El_injection_rates->setrate(Event_map,El_injection_events[Event_map]->rate);
        }
        else {
            Event_map = inject_node;
        }
        
        Ho_injection_events[Event_map]->Set_injection_event(right_electrode, inject_node, Hole, formalism, 0.0, lrto, globevent);
        Ho_injection_rates->setrate(Event_map,Ho_injection_events[Event_map]->rate);
    }
}

void Events::Initialize_eventvector_for_device(bool dual_injection, vector <Carrier*> electrons, vector <Carrier*> holes, Node* left_electrode, Node* right_electrode, int maxpairdegree){ //
    
    El_non_injection_events.clear();
    Ho_non_injection_events.clear();
    Grow_non_injection_eventvector(electrons.size(), electrons, El_non_injection_events, maxpairdegree);
    Grow_non_injection_eventvector(holes.size(), holes,Ho_non_injection_events, maxpairdegree);

    
    if(dual_injection) {
        El_injection_events.clear();
        Ho_injection_events.clear();
        Initialize_injection_eventvector(left_electrode,El_injection_events);
        Initialize_injection_eventvector(right_electrode,El_injection_events);
        Initialize_injection_eventvector(left_electrode,Ho_injection_events);
        Initialize_injection_eventvector(right_electrode,Ho_injection_events);
    }
    else {
        El_injection_events.clear();
        Ho_injection_events.clear();
        Initialize_injection_eventvector(left_electrode,El_injection_events);
        Initialize_injection_eventvector(right_electrode,Ho_injection_events); 
    }
    
}

void Events::Initialize_injection_eventvector(Node* electrode, vector<Event*> eventvector){

    for (int inject_node = 0; inject_node<electrode->pairing_nodes.size(); inject_node++) {

        Event *newEvent = new Event();
        eventvector.push_back(newEvent);
        
    } 
}

void Events::Grow_non_injection_eventvector(int carrier_grow_size, vector<Carrier*> carriers, vector<Event*> eventvector,int maxpairdegree){
    
    int old_nr_carriers = div(eventvector.size(),maxpairdegree).quot; //what was the number of carriers that we started with?
    
    for(int carrier_ID = old_nr_carriers; carrier_ID<old_nr_carriers+carrier_grow_size; carrier_ID++) {
        for(int jump_ID = 0; jump_ID<maxpairdegree;jump_ID++) {

            Event *newEvent = new Event();
            eventvector.push_back(newEvent);

            newEvent->carrier = carriers[carrier_ID];
        }         
    }    
}


}} 

#endif
