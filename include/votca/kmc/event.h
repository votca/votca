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

#ifndef __VOTCA_KMC_EVENT_H_
#define __VOTCA_KMC_EVENT_H_

#include <votca/kmc/state.h>
#include <votca/kmc/graph.h>
#include <votca/kmc/globaleventinfo.h>

namespace votca { namespace kmc {
  
using namespace std;

enum From_step_event {Fromtransfer, Injection, Fromnotinbox};
enum To_step_event {Totransfer, Recombination, Collection, Blocking, Tonotinbox};

/*
 * Abstract base class for all events  
 */
class Event {
    
public:

    From_step_event fromtype;
    To_step_event totype;
    double rate;
    Graph* graph;
    Carrier* carrier;
    Globaleventinfo* globaleventinfo;
    
    void Set_injection_event(Node* electrode, int injectnode_ID, CarrierType carrier_type, string formalism,
                              double from_longrange, double to_longrange, Globaleventinfo* globaleventinfo);
    void Set_non_injection_event(vector<Node*> nodes, Carrier* carrier, int jump_ID, string formalism, 
                                 double from_longrange, double to_longrange, Globaleventinfo* globaleventinfo);
    
    From_step_event Determine_non_injection_from_event_type(Carrier* carrier);
    To_step_event Determine_non_injection_to_event_type(Carrier* carrier, int jumpID, Node* carriernode);
    To_step_event Determine_injection_to_event_type(CarrierType carrier_type, Node* electrode, int inject_nodeID);

    double Compute_event_rate(Node* fromnode, int jump_ID, CarrierType carrier_type,
                            From_step_event from_event_type, To_step_event to_event_type, string formalism,
                            double from_shortrange, double to_shortrange, double from_longrange, double to_longrange,
                            Globaleventinfo* globaleventinfo);
    
};

void Event::Set_injection_event(Node* electrode, int injectnode_ID, CarrierType carrier_type, string formalism,
                              double from_longrange, double to_longrange, Globaleventinfo* globaleventinfo) {
    
    fromtype = Injection;
    totype = Determine_injection_to_event_type(carrier_type, electrode, injectnode_ID);
    rate = Compute_event_rate(electrode, injectnode_ID, carrier_type, fromtype, totype, formalism,
                              0, 0.0, from_longrange, to_longrange, globaleventinfo);
       
}

void Event::Set_non_injection_event(vector<Node*> nodes, Carrier* carrier, int jump_ID, string formalism, 
                                 double from_longrange, double to_longrange, Globaleventinfo* globaleventinfo) {
    
    fromtype = Determine_non_injection_from_event_type(carrier);
    totype = Determine_non_injection_to_event_type(carrier, jump_ID, nodes[carrier->carrier_node_ID]);
    rate = Compute_event_rate(nodes[carrier->carrier_node_ID], jump_ID, carrier->carrier_type, fromtype, totype, formalism,
                                     carrier->srfrom, carrier->srto[jump_ID], from_longrange, to_longrange,
                                     globaleventinfo);    
}

double Event::Compute_event_rate(Node* fromnode, int jump_ID, CarrierType carrier_type,
                                     From_step_event from_event_type, To_step_event to_event_type, string formalism,
                                     double from_shortrange, double to_shortrange, double from_longrange, double to_longrange,
                                     Globaleventinfo* globevent){

    Node* jumptonode = fromnode->pairing_nodes[jump_ID];

    double prefactor = 1.0;
    double charge;
    double static_node_energy_from;
    double static_node_energy_to;
    
    if(carrier_type == Electron) {
        charge = 1.0;
        prefactor *= globevent->electron_prefactor;
        static_node_energy_from = fromnode->static_electron_node_energy;
        static_node_energy_to = jumptonode->static_electron_node_energy; 
    }
    else {
        charge = -1.0;
        prefactor *= globevent->hole_prefactor;
        static_node_energy_from = fromnode->static_hole_node_energy;
        static_node_energy_to = jumptonode->static_hole_node_energy;
    }
    
    //first calculate quantum mechanical wavefunction overlap
    myvec distancevector = fromnode->static_event_info[jump_ID].distance;
    double distance = abs(distancevector);
    
    double distancefactor = exp(-1.0*globevent->alpha*distance);
  
    //second, calculate boltzmann factor (Coulomb interaction still have to be implemented)
   
    double init_energy;
    double final_energy;
    double selfimpot_from = fromnode->self_image_potential;
    double selfimpot_to = jumptonode->self_image_potential;
    double fromtype_energy = 0.0;
    double totype_energy = 0.0;
    
    if(from_event_type == Injection) {
        fromtype_energy -= globevent->injection_barrier;
        prefactor *= globevent->injection_prefactor;
    }
    if(to_event_type == Recombination) {
        totype_energy -= globevent->binding_energy;
        prefactor *= globevent->recombination_prefactor;
    }
    if(to_event_type == Collection) {
        totype_energy -= globevent->injection_barrier;
        prefactor *= globevent->collection_prefactor;
    }
    
    double coulomb_from;
    double coulomb_to;
    
    if ((from_event_type != Injection)) { 
        coulomb_from = globevent->coulomb_strength*(from_shortrange+charge*from_longrange);
        coulomb_to = globevent->coulomb_strength*(to_shortrange+charge*to_longrange);
    }
    else if(from_event_type == Injection){
        coulomb_from = 0.0;
        coulomb_to = charge*globevent->coulomb_strength*(jumptonode->injection_potential+to_longrange);
    }
      
    init_energy = static_node_energy_from + selfimpot_from + coulomb_from;
    final_energy = static_node_energy_to + selfimpot_to + coulomb_to;
    
    double energycontrib;
    double energyfactor;
    
    if (formalism == "Miller") {
        if(to_event_type == Blocking) {
            energyfactor = 0.0; // Keep this here for eventual simulation of bipolaron formation for example
        }
        else if((from_event_type == Fromnotinbox)&&(to_event_type == Tonotinbox)) {
            energyfactor = 0.0; // Keep this here for eventual simulation of one-site events (on-node generation)
        }
        else {
            energycontrib = final_energy - init_energy -charge*globevent->efield*distancevector.x();
            if (energycontrib>0.0) {
                energyfactor = exp(-1.0*globevent->beta*energycontrib);
            }
            else {
                energyfactor = 1.0;
            }
        }
    }
    
    double jump_rate = prefactor*distancefactor*energyfactor;
    return jump_rate;    
}

From_step_event Event::Determine_non_injection_from_event_type(Carrier* carrier){
    
    From_step_event from_type;
    
    if(!carrier->is_in_sim_box){
        from_type = Fromnotinbox;
    }
    else {
        from_type = Fromtransfer;
    }
    
    return from_type;
    
}

To_step_event Event::Determine_non_injection_to_event_type(Carrier* carrier, int jumpID, Node* carriernode){
    
    To_step_event to_type;
    
    if(!carrier->is_in_sim_box){
        to_type = Tonotinbox;
    }
    else {
        int node_degree = carriernode->pairing_nodes.size();
        if(jumpID < node_degree) { // hopping event exists in graph
            Node* jumpnode = carriernode->pairing_nodes[jumpID];
            if(jumpnode->carriers_on_node.empty()){
                to_type = Totransfer;
            }
            else if(jumpnode->carriers_on_node[0]->carrier_type == carrier->carrier_type) {
                to_type = Blocking;
            }
            else if(jumpnode->carriers_on_node[0]->carrier_type != carrier->carrier_type) {
                to_type = Recombination;
            }
        }
        else {
            to_type = Tonotinbox;
        }
    }
    
    return to_type;
}

To_step_event Event::Determine_injection_to_event_type(CarrierType carrier_type, Node* electrode, int inject_nodeID){
    
    Node* injectnode = electrode->pairing_nodes[inject_nodeID];
    if(injectnode->carriers_on_node.empty()){
        totype = Totransfer;
    }
    else if(injectnode->carriers_on_node[0]->carrier_type == carrier_type) {
        totype = Blocking;
    }
    else if(injectnode->carriers_on_node[0]->carrier_type != carrier_type) {
        totype = Recombination;
    }
    
    return totype;
}



}} 

#endif

