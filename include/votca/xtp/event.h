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

#ifndef __VOTCA_KMC_EVENT_H_
#define __VOTCA_KMC_EVENT_H_

#include <votca/xtp/carrier.h>
#include <votca/xtp/graphkmc.h>
#include <votca/xtp/eventinfo.h>
#include <votca/xtp/state.h>
#include <votca/xtp/longrange.h>
#include <votca/xtp/statereservoir.h>

namespace votca { namespace xtp {
  
using namespace std;

enum Final_Event_Type {TransferTo, Collection, CollectiontoNode, Recombination, Blocking, Notinbox};
enum Init_Event_Type {Injection, InjectionfromNode, TransferFrom, Notinboxfrom};
enum Action{Add, Remove, None };
enum Action_pair{Transfer, Other };

class Event {
    
public:
    
    Event(int id, int type) {
        _id = id;
        _init_type = (int) Notinboxfrom;
        _final_type = type;
        _rate = 0.0;
    }

    Event(int id, Link* link, int carrier_type, StateReservoir* state, Longrange* longrange, Eventinfo* eventinfo){
        _id = id;
        if(link->node1()->type() != (int) NormalNode ) _injection_potential = 0.0;
        Set_event(link, carrier_type, state, longrange, eventinfo);
    }
    
    double &rate() { return _rate; }
    int &init_type() { return _init_type; }
    int &final_type() { return _final_type;}
    Link* link() {return _link;}
    int &id() { return _id;}
    int &carrier_type() { return _carrier_type;}
    int &action_node1() { return _action_node1;}
    int &action_node2() { return _action_node2;}
    int &action_pair()  { return _action_pair; }
    double &transferfactor() { return _transferfactor;}
    double &energyfactor() { return _energyfactor;}
    
    double self_impot1(Eventinfo* eventinfo) { return Determine_self_coulomb(_link->node1(), eventinfo);}
    double self_impot2(Eventinfo* eventinfo) { return Determine_self_coulomb(_link->node2(), eventinfo);}
    double lr_coulomb1(Longrange* longrange, Eventinfo* eventinfo) { return Determine_lr_coulomb(_link->node1(), longrange, eventinfo);}
    double lr_coulomb2(Longrange* longrange, Eventinfo* eventinfo) { return Determine_lr_coulomb(_link->node2(), longrange, eventinfo);}
    double sr_coulomb1(StateReservoir* state, Eventinfo* eventinfo) { return Determine_from_sr_coulomb(_link->node1(), state, eventinfo);}
    double sr_coulomb2(StateReservoir* state, Eventinfo* eventinfo) { return Determine_to_sr_coulomb(_link->node1(), state, eventinfo);}
    
    void Set_event(Link* link, int carrier_type, StateReservoir* state, Longrange* longrange, Eventinfo* eventinfo);

    /// Determine rate
    void Determine_rate(StateReservoir* state, Longrange* longrange, Eventinfo* eventinfo);
    /// Set rate to value
    void Set_rate(double rate) {_rate = rate;}
    /// Set out of box
    void Set_not_in_box_event() {_final_type = Notinbox; _rate = 0.0;}
    /// Set injection potential
    void Set_injection_potential(double injection_potential) {_injection_potential = injection_potential;}
    /// Add injection potential
    void Add_injection_potential(double injection_potential) {_injection_potential += injection_potential;}
    
    /// Determine initial event type
    int Determine_init_event_type(Node* node1);
    /// Determine non injection event type in case two carriers are on linking nodes
    int Determine_final_event_type(int carrier_type1, int carrier_type2, Node* node1, Node* node2, Eventinfo* eventinfo);
    /// Determine non injection event type in case only one carrier is on a link
    int Determine_final_event_type(Node* node1, Node* node2);
    /// Determine action flag for node 1
    int Determine_action_flag_node1();
    /// Determine action flag for node 2
    int Determine_action_flag_node2();
    
    int Determine_action_flag_pair();

    /// Determine short range coulomb potential at node from which hop occur
    double Determine_from_sr_coulomb(Node* node, StateReservoir* state, Eventinfo* eventinfo); 
    
    /// Determine short range coulomb potential at node to which hop occur
    double Determine_to_sr_coulomb(Node* node, StateReservoir* state, Eventinfo* eventinfo);
    
    /// Determine long range coulomb potential for given node (equal to 0 in case of electrode node)
    double Determine_lr_coulomb(Node* node, Longrange* longrange, Eventinfo* eventinfo);
 
    double Determine_self_coulomb(Node* node, Eventinfo* eventinfo);    
    
private:
                         
    Link* _link;
    int _carrier_type;
    
    double _rate;
    double _energyfactor;
    double _transferfactor;
    
    int _final_type;
    int _init_type;
    int _id;
    
    int _action_node1;
    int _action_node2;
    int _action_pair;
    
    //int _layer_node1;
    //int _layer_node2;
    
    double _injection_potential;
    
    
};



}} 

#endif

