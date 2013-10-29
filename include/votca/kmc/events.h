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
#include <votca/kmc/longrange.h>
#include <votca/kmc/bsumtree.h>

namespace votca { namespace kmc {
  
using namespace std;

/*
 * Abstract base class for all events  
 */
class Events {
    
public:

    vector<Event*> El_non_injection_events;
    vector<Event*> Ho_non_injection_events;
    vector<Event*> El_injection_events;
    vector<Event*> Ho_injection_events;
    State* state;
    Graph* graph;

    void Set_non_injection_event(Carrier* carrier, int jumpID);
    void Set_injection_event(Node* electrode, int inject_nodeID, CarrierType carrier_type, bool dual_injection);

    From_step_event Determine_non_injection_from_event_type(Carrier* carrier);
    To_step_event Determine_non_injection_to_event_type(Carrier* carrier, int jumpID);
    To_step_event Determine_injection_to_event_type(CarrierType carrier_type, Node* electrode, int inject_nodeID);
    
    double Compute_non_injection_event_rate(Carrier* carrier, int jumpID, From_step_event from, To_step_event to);
    double Compute_injection_event_rate(CarrierType carrier_type, Node* electrode, int inject_nodeID, From_step_event from, To_step_event to);  
    
    void Recompute_all_injection_events(CarrierType carrier_type, bool dual_injection);
    void Recompute_all_non_injection_events(CarrierType carrier_type);
    
    void Initialize_events_for_device(bool dual_injection); // if dual_injection is true, both types of carriers are injected from both electrodes
    void Grow_non_injection_events(int event_grow_size, CarrierType carrier_type);
    
private:
    void Initialize_injection_events(Node* electrode, CarrierType carrier_type);
};

/*double Events::compute_non_injection_event_rate(Carrier* carrier, int jumpID, From_step_event from, To_step_event to) {
    
}*/

void Events::Set_non_injection_event(Carrier* carrier, int jumpID) {
    
    Event* non_inject_event;
    int Event_map = carrier->carrier_ID*graph->max_pair_degree + jumpID;
    
    if(carrier->carrier_type == Electron) {
        non_inject_event = El_non_injection_events[Event_map];
    }
    else if(carrier->carrier_type == Hole) {
        non_inject_event = Ho_non_injection_events[Event_map];
    }
    
    From_step_event fromtype = Determine_non_injection_from_event_type(carrier);
    To_step_event totype = Determine_non_injection_to_event_type(carrier, jumpID);
    
    non_inject_event->fromtype = fromtype;
    non_inject_event->totype = totype;
    non_inject_event->rate = Compute_non_injection_event_rate(carrier, jumpID, fromtype, totype);    
    
}

void Events::Set_injection_event(Node* electrode, int inject_nodeID, CarrierType carrier_type, bool dual_injection) {
    
    Event* inject_event;
    int Event_map;
    int electrode_ID;
    
    if(electrode->node_type = LeftElectrode) {electrode_ID = 0;}
    if(electrode->node_type = RightElectrode) {electrode_ID = 1;}
    
    if(dual_injection) {
        Event_map = electrode_ID*graph->left_electrode->pairing_nodes.size() + inject_nodeID;
    }
    else {
        Event_map = inject_nodeID;
    }
    
    if(carrier_type == Electron) {
        inject_event = El_injection_events[Event_map];
    }
    else if(carrier_type == Hole) {
        inject_event = Ho_injection_events[Event_map];
    }

    From_step_event fromtype = Injection;
    To_step_event totype = Determine_injection_to_event_type(carrier_type, electrode, inject_nodeID);
    
    inject_event->fromtype = fromtype;
    inject_event->totype = totype;
    inject_event->rate = Compute_injection_event_rate(carrier_type, electrode, inject_nodeID, fromtype, totype);    
    
}

From_step_event Events::Determine_non_injection_from_event_type(Carrier* carrier){
    
    From_step_event from;
    if(!carrier->is_in_sim_box){
        from = Fromnotinbox;
    }
    else {
        from = Fromtransfer;
    }
    
    return from;
}

To_step_event Events::Determine_non_injection_to_event_type(Carrier* carrier, int jumpID){
    
    To_step_event to;
    if(!carrier->is_in_sim_box){
        to = Tonotinbox;
    }
    else {
        if(jumpID<graph->nodes[carrier->carrier_node_ID]->pairing_nodes.size()) { // hopping event exists in graph
            Node* jumpnode = graph->nodes[carrier->carrier_node_ID]->pairing_nodes[jumpID];
            if(jumpnode->carriers_on_node.empty()){
                to = Totransfer;
            }
            else if(jumpnode->carriers_on_node[0]->carrier_type == carrier->carrier_type) {
                to = Blocking;
            }
            else if(jumpnode->carriers_on_node[0]->carrier_type != carrier->carrier_type) {
                to = Recombination;
            }
        }
        else {
            to = Tonotinbox;
        }
    }
    
    return to;
}

To_step_event Events::Determine_injection_to_event_type(CarrierType carrier_type, Node* electrode, int inject_nodeID){
    
    To_step_event to;
    Node* injectnode = electrode->pairing_nodes[inject_nodeID];
    if(injectnode->carriers_on_node.empty()){
        to = Totransfer;
    }
    else if(injectnode->carriers_on_node[0]->carrier_type == carrier_type) {
        to = Blocking;
    }
    else if(injectnode->carriers_on_node[0]->carrier_type != carrier_type) {
        to = Recombination;
    }

    return to;
}

void Events::Recompute_all_non_injection_events(CarrierType carrier_type) {
    
    vector<Carrier*> carriers;
    if(carrier_type == Electron) {
        carriers = state->electrons;
    }
    else if(carrier_type == Hole) {
        carriers = state->holes;
    }
    
    for (int carrier_ID = 0; carrier_ID<carriers.size(); carrier_ID++) {
        for (int ipair = 0; ipair < graph->max_pair_degree;ipair++){
            Set_non_injection_event(carriers[carrier_ID], ipair);
        }
    }
}

void Events::Recompute_all_injection_events(CarrierType carrier_type, bool dual_injection) {
    for (int inject_node = 0; inject_node<graph->left_electrode->pairing_nodes.size(); inject_node++) {
        Set_injection_event(graph->left_electrode, inject_node, carrier_type, dual_injection);
    }
    
    for (int inject_node = 0; inject_node<graph->right_electrode->pairing_nodes.size(); inject_node++) {
        Set_injection_event(graph->right_electrode, inject_node, carrier_type, dual_injection);
    }
}

void Events::Initialize_events_for_device(bool dual_injection){ //
    
    Grow_non_injection_events(state->electrons.size(),Electron);
    Grow_non_injection_events(state->holes.size(),Hole);
    
    if(dual_injection) {
        El_injection_events.clear();
        Ho_injection_events.clear();    
        Initialize_injection_events(graph->left_electrode,Electron);
        Initialize_injection_events(graph->right_electrode,Electron);
        Initialize_injection_events(graph->left_electrode,Hole);
        Initialize_injection_events(graph->right_electrode,Hole);        
    }
    else {
        El_injection_events.clear();
        Ho_injection_events.clear();  
        Initialize_injection_events(graph->left_electrode,Electron);
        Initialize_injection_events(graph->right_electrode,Hole); 
    }
}

void Events::Initialize_injection_events(Node* electrode, CarrierType carrier_type){
    for (int inject_node = 0; inject_node<electrode->pairing_nodes.size(); inject_node++) {

        vector<Carrier*> carriers;
        vector<Event*> events;
        
        if(carrier_type == Electron) {
            carriers = state->electrons;
            events = El_injection_events;
        }
        else if(carrier_type == Hole) {
            carriers = state->holes;
            events = Ho_injection_events;
        }
        
        Event *newEvent = new Event();
        events.push_back(newEvent);

        From_step_event fromtype = Injection;
        To_step_event totype = Determine_injection_to_event_type(carrier_type, electrode, inject_node); 
        
        newEvent->fromtype = fromtype;
        newEvent->totype = totype;        
        newEvent->rate = Compute_injection_event_rate(carrier_type, electrode, inject_node, fromtype, totype);
    } 
}    

void Events::Grow_non_injection_events(int carrier_grow_size, CarrierType carrier_type){

    vector<Carrier*> carriers;
    vector<Event*> events;
    
    if(carrier_type == Electron) {
        carriers = state->electrons;
        events = El_non_injection_events;
    }
    else if(carrier_type == Hole) {
        carriers = state->holes;
        events = Ho_non_injection_events;
    }
    
    int old_nr_carriers = div(events.size(),graph->max_pair_degree).quot;
    
    for(int carrier_ID = old_nr_carriers; carrier_ID<old_nr_carriers+carrier_grow_size; carrier_ID++) {
        for(int jump_ID = 0; jump_ID<graph->max_pair_degree;jump_ID++) {
            
            Event *newEvent = new Event();
            events.push_back(newEvent);

            From_step_event fromtype = Determine_non_injection_from_event_type(carriers[carrier_ID]);
            To_step_event totype = Determine_non_injection_to_event_type(carriers[carrier_ID], jump_ID);            
            
            newEvent->fromtype = fromtype;
            newEvent->totype = totype; 
            newEvent->rate = Compute_non_injection_event_rate(carriers[carrier_ID], jump_ID, fromtype, totype);
        }         
    }    
}


}} 

#endif

