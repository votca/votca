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

#ifndef __VOTCA_KMC_OHMIC_EVENT_H_
#define __VOTCA_KMC_OHMIC_EVENT_H_

#include <votca/kmc/carrier.h>
#include <votca/kmc/graphdevice.h>
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/state.h>
#include <votca/kmc/longrange.h>

namespace votca { namespace kmc {
  
using namespace std;

class Ohmic_event {
    
public:
    
    Ohmic_event(int id, Node* node, int electrode, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
        _id = id;
        Set_ohmic_event(node, electrode, state, longrange, eventinfo);
    }
 
    double &prob() { return _prob; }
    Node* node() {return _node;}
    int &id() { return _id;}
    
    void Set_ohmic_event(Node* node, int electrode, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    
    void Determine_prob(int electrode, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    
    /// Set injection potential
    void Set_injection_potential(double injection_potential) {_injection_potential = injection_potential;}
    /// Add injection potential
    void Add_injection_potential(double injection_potential) {_injection_potential += injection_potential;}
   
    
protected:
                         
    Node* _node;
    int _carrier_type;
    
    double _prob;
    int _id;
    
    double _injection_potential;
    
};

void Ohmic_event::Determine_prob(int electrode, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
    
    double charge;
    double static_node_energy;
    double static_electrode_energy = eventinfo->avholeenergy;
    const double kB   = 8.617332478E-5; // ev/K     
    
    // we are interested in the difference with the electrode

    if(_carrier_type == (int) Electron) {
        charge = -1.0;
        static_node_energy = dynamic_cast<NodeSQL*>(_node)->eCation() + dynamic_cast<NodeSQL*>(_node)->UcCnNe();
    }
    else if(_carrier_type == (int) Hole) {
        charge = 1.0;
        static_node_energy = dynamic_cast<NodeSQL*>(_node)->eAnion() + dynamic_cast<NodeSQL*>(_node)->UcCnNh();
    }    
    
    double electrode_energy;
    double node_energy;
    
    double inject_bar = 1.0*eventinfo->injection_barrier;
    double sr_coulomb = eventinfo->coulomb_strength*_injection_potential;
    double lr_coulomb = eventinfo->coulomb_strength*charge*longrange->Get_cached_longrange(dynamic_cast<NodeDevice*>(_node)->layer());
    double selfimpot = dynamic_cast<NodeDevice*>(_node)->self_image();    

    electrode_energy = static_electrode_energy;
    node_energy = static_node_energy + selfimpot + inject_bar + sr_coulomb + lr_coulomb;
 
    double distance_to_electrode;
    votca::tools::vec node_pos = _node->position();
    
    if(electrode==0) { distance_to_electrode = node_pos.x(); }
    else { distance_to_electrode = eventinfo->simboxsize.x() - node_pos.x(); }
 
    double bar = node_energy - electrode_energy -charge*(eventinfo->efield_x*distance_to_electrode);
    
    _prob = 1.0/(1.0+exp(bar/(kB*eventinfo->temperature)));
}

void Ohmic_event::Set_ohmic_event(Node* node, int electrode, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
    
    _node = node;
    _carrier_type = eventinfo->ohmic_cartype;
    
    Determine_prob(electrode, state, longrange, eventinfo);

}

}} 

#endif

