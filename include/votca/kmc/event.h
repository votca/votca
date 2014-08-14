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

#include <votca/kmc/carrier.h>
#include <votca/kmc/graphkmc.h>
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/state.h>
#include <votca/kmc/longrange.h>

namespace votca { namespace kmc {
  
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
    inline int Determine_init_event_type(Node* node1);
    /// Determine non injection event type in case two carriers are on linking nodes
    inline int Determine_final_event_type(int carrier_type1, int carrier_type2, Node* node1, Node* node2, Eventinfo* eventinfo);
    /// Determine non injection event type in case only one carrier is on a link
    inline int Determine_final_event_type(Node* node1, Node* node2);
    /// Determine action flag for node 1
    inline int Determine_action_flag_node1();
    /// Determine action flag for node 2
    inline int Determine_action_flag_node2();
    
    inline int Determine_action_flag_pair();

    /// Determine short range coulomb potential at node from which hop occur
    inline double Determine_from_sr_coulomb(Node* node, StateReservoir* state, Eventinfo* eventinfo); 
    
    /// Determine short range coulomb potential at node to which hop occur
    inline double Determine_to_sr_coulomb(Node* node, StateReservoir* state, Eventinfo* eventinfo);
    
    /// Determine long range coulomb potential for given node (equal to 0 in case of electrode node)
    inline double Determine_lr_coulomb(Node* node, Longrange* longrange, Eventinfo* eventinfo);
 
    inline double Determine_self_coulomb(Node* node, Eventinfo* eventinfo);    
    
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
    
    int _layer_node1;
    int _layer_node2;
    
    double _injection_potential;
    
    
};

inline int Event::Determine_final_event_type(Node* node1, Node* node2) {
    if (node2->type() == (int) NormalNode){                                                                  return (int) TransferTo;}
    else if (node2->type() == (int) LeftElectrodeNode || node2->type() == (int) RightElectrodeNode){         return (int) Collection;} // Collection at electrode
}

inline int Event::Determine_final_event_type(int carrier_type1, int carrier_type2, Node* node1, Node* node2, Eventinfo* eventinfo) {
    if (((carrier_type1 == (int) Electron) && (carrier_type2 == (int) Electron)) || ((carrier_type1 == (int) Hole) && (carrier_type2 == (int) Hole))){
        if(!eventinfo->no_blocking) {return (int) Blocking;} else { return (int) TransferTo;} // Blocking
    }
    else if (((carrier_type1 == (int) Electron) && (carrier_type2 == (int) Hole))|| ((carrier_type1 == (int) Hole) && (carrier_type2 == (int) Electron))){ 
        return (int) Recombination; // Recombination
    }
}

inline int Event::Determine_init_event_type(Node* node1) {
    if(node1->type() == (int) NormalNode){                                                                   return (int) TransferFrom;}
    else if((node1->type() == (int) LeftElectrodeNode) || (node1->type() == (int) RightElectrodeNode)){      return (int) Injection;   }
}

inline int Event::Determine_action_flag_pair() {
    int action_pair;
    if(_init_type == TransferFrom && _final_type == TransferTo) {action_pair = (int) Transfer; }
    else                                                        {action_pair = (int) Other;    }
    return action_pair;
}

inline int Event::Determine_action_flag_node1() {
    int action_node1;
    if(_init_type == Injection)           {action_node1 = (int) None;   } // injection
    else if(_init_type == TransferFrom)   {action_node1 = (int) Remove; } // transfer
    return action_node1;    
}

inline int Event::Determine_action_flag_node2() {
    int action_node2;
    Node* node2 = _link->node2();
    if(_final_type == TransferTo)            {action_node2 = (int) Add;    } // transfer
    else if(_final_type == Collection)       {action_node2 = (int) None;   } // collection
    else if(_final_type == Recombination)    {action_node2 = (int) Remove; } // recombination
    return action_node2;    
}

inline double Event::Determine_from_sr_coulomb(Node* node, StateReservoir* state, Eventinfo* eventinfo) {

    double coulomb_from;
    if(_init_type == Injection)   {coulomb_from = 0.0;                                                                           }
    else                          {coulomb_from = eventinfo->coulomb_strength*state->GetCarrier(node->occ())->on_site_coulomb(); }
    return coulomb_from;
}

inline double Event::Determine_to_sr_coulomb(Node* node, StateReservoir* state, Eventinfo* eventinfo) {

    double coulomb_to;
    if(_final_type == Collection) {coulomb_to = 0.0;                                                                                       }
    else if(_init_type != Injection) {coulomb_to = eventinfo->coulomb_strength*state->GetCarrier(node->occ())->to_site_coulomb(_link->id());}
    else if(_init_type == Injection) { coulomb_to = eventinfo->coulomb_strength*_injection_potential; }
    return coulomb_to;
}

inline double Event::Determine_lr_coulomb(Node* node, Longrange* longrange, Eventinfo* eventinfo) {
    double lr_coulomb;

    if(node->type() == (int) NormalNode) {
        if(!eventinfo->longrange_slab) {
            lr_coulomb = eventinfo->lr_coulomb_strength*longrange->Get_cached_longrange(dynamic_cast<NodeDevice*>(node)->layer());
        }
        else {
            lr_coulomb = eventinfo->lr_coulomb_strength*longrange->Get_cached_longrange_slab(node->id());
        }
    }
    else {lr_coulomb = 0.0;} //potential at electrodes = 0.0;
    return lr_coulomb;
}

inline double Event::Determine_self_coulomb(Node* node, Eventinfo* eventinfo){
        double self_coulomb;
        if(node->type() == (int) NormalNode) {self_coulomb = eventinfo->coulomb_strength*eventinfo->self_image_prefactor*dynamic_cast<NodeDevice*>(node)->self_image();}
        else {self_coulomb = 0.0;}
        return self_coulomb;
}

void Event::Determine_rate(StateReservoir* state, Longrange* longrange, Eventinfo* eventinfo) {

    Node* node1 = _link->node1();
    votca::tools::vec node1vec = node1->position();
    double leftnode1pos = node1vec.z();
    double rightnode1pos = eventinfo->simboxsize.z() - leftnode1pos;
    Node* node2 = _link->node2();
    votca::tools::vec node2vec = node2->position();
    double leftnode2pos = node2vec.z();
    double rightnode2pos = eventinfo->simboxsize.z() - leftnode2pos;
    
    double prefactor = 1.0; // total prefactor
    double charge;
    double static_node_energy_from;
    double static_node_energy_to;

    const double hbar = 6.58211928E-16; // eV*s
    const double Pi   = 3.14159265358979323846264338327950288419716939937510;
    const double kB   = 8.617332478E-5; // ev/K     

    if(_carrier_type == (int) Electron) {
        charge = -1.0;
        prefactor = prefactor*(eventinfo->electron_transport_prefactor);
        static_node_energy_from = dynamic_cast<NodeSQL*>(node1)->eAnion() + dynamic_cast<NodeSQL*>(node1)->UcCnNe();
        static_node_energy_to = dynamic_cast<NodeSQL*>(node2)->eAnion() + dynamic_cast<NodeSQL*>(node2)->UcCnNe();
        if(_init_type == Injection) static_node_energy_from = eventinfo->avelectronenergy; 
        if(_final_type == Collection) static_node_energy_to = eventinfo->avelectronenergy;
    }
    else if(_carrier_type == (int) Hole) {
        charge = 1.0;
        double conversion1 = 0.57;
        double conversion2 = 0.57;
        //if(eventinfo->novikov) {
        //    conversion1 = (1.0 - (0.9/(2.0*leftnode1pos))*(1.0 - exp(-2.0*leftnode1pos/0.9)))*(1.0 - (0.9/(2.0*rightnode1pos))*(1.0 - exp(-2.0*rightnode1pos/0.9)));
        //    conversion2 = (1.0 - (0.9/(2.0*leftnode2pos))*(1.0 - exp(-2.0*leftnode2pos/0.9)))*(1.0 - (0.9/(2.0*rightnode2pos))*(1.0 - exp(-2.0*rightnode2pos/0.9)));
        //}
        prefactor = prefactor*(eventinfo->hole_transport_prefactor);
        double temp_static_node_energy_from = dynamic_cast<NodeSQL*>(node1)->eCation() + dynamic_cast<NodeSQL*>(node1)->UcCnNh();
        static_node_energy_from = temp_static_node_energy_from;
        double temp_static_node_energy_to = dynamic_cast<NodeSQL*>(node2)->eCation() + dynamic_cast<NodeSQL*>(node2)->UcCnNh();
        static_node_energy_to = temp_static_node_energy_to;
        if(_init_type == Injection) {
            static_node_energy_from = eventinfo->avholeenergy; 
            static_node_energy_to = conversion2*(temp_static_node_energy_to - eventinfo->avholeenergy) + eventinfo->avholeenergy;
        }
        if(_final_type == Collection) {
            static_node_energy_from = conversion1*(temp_static_node_energy_from - eventinfo->avholeenergy) + eventinfo->avholeenergy;
            static_node_energy_to = eventinfo->avholeenergy;
        }
    }

    //first transfer integrals

    double _transferfactor = 1.0;
    votca::tools::vec distancevector = _link->r12();
    double Reorg;
    double Jeff2;
    double distance;

    if (eventinfo->formalism == "Miller") {
        distance = abs(distancevector);
        _transferfactor = exp(-2.0*eventinfo->alpha*distance);
    }
    else if((eventinfo->formalism == "Marcus")&&(_init_type==Injection||_final_type==Collection)) {
        distance = abs(distancevector);
        Jeff2 = exp(-2.0*eventinfo->alpha*distance);

        // take Reorg equal to inside organic material (obviously not true))
        if(_carrier_type == (int) Electron) {
            Reorg = eventinfo->electron_injection_reorg;
        }
        if(_carrier_type == (int) Hole) {
            Reorg = eventinfo->hole_injection_reorg;
        }
        _transferfactor = (2*Pi/hbar)*(Jeff2/sqrt(4*Pi*Reorg*kB*eventinfo->temperature));
    }
    else if(eventinfo->formalism == "Marcus") {
        if(_carrier_type == (int) Electron) {
            Jeff2 = dynamic_cast<LinkSQL*>(_link)->Jeff2e();
            Reorg = dynamic_cast<NodeSQL*>(node1)->UnCnNe() + dynamic_cast<NodeSQL*>(node2)->UcNcCe() + dynamic_cast<LinkSQL*>(_link)->lOe();
        }
        if(_carrier_type == (int) Hole) {
            Jeff2 = dynamic_cast<LinkSQL*>(_link)->Jeff2h();
            Reorg = dynamic_cast<NodeSQL*>(node1)->UnCnNh() + dynamic_cast<NodeSQL*>(node2)->UcNcCh() + dynamic_cast<LinkSQL*>(_link)->lOh();
        }
        _transferfactor = (2*Pi/hbar)*(Jeff2/sqrt(4*Pi*Reorg*kB*eventinfo->temperature));
    }

    //second, calculate boltzmann factor

    double init_energy;
    double final_energy;
    double from_event_energy = 0.0;
    double to_event_energy = 0.0;

    if(_init_type == Injection)           { 
        if(node1->type() == (int) LeftElectrodeNode) {from_event_energy -= eventinfo->left_injection_barrier;}
        else if(node1->type() == (int) RightElectrodeNode) {from_event_energy -= eventinfo->right_injection_barrier;}
        prefactor *= eventinfo->injection_prefactor;    
    } // injection
    else if(_init_type == TransferFrom)   {                                                                                                    } // transfer

    if(_final_type == TransferTo)         {                                                                                                    } // transfer
    else if(_final_type == Collection)    { 
        if(node2->type() == (int) LeftElectrodeNode) {to_event_energy -= eventinfo->left_injection_barrier;}
        else if(node2->type() == (int) RightElectrodeNode) {to_event_energy -= eventinfo->right_injection_barrier;}
        prefactor *= eventinfo->collection_prefactor;               
    } // collection
    else if(_final_type == Recombination) { to_event_energy   -= eventinfo->binding_energy;    prefactor *= eventinfo->recombination_prefactor;} // recombination

    double sr_coulomb_from;
    double sr_coulomb_to;
    if(!eventinfo->norc) {
        sr_coulomb_from = Determine_from_sr_coulomb(node1, state, eventinfo);
        sr_coulomb_to = Determine_to_sr_coulomb(node1, state, eventinfo);
    }
    else {
        sr_coulomb_from = 0.0;
        sr_coulomb_to = 0.0;
    }
    double selfimpot_from;
    double selfimpot_to;
    double lr_coulomb_from;
    double lr_coulomb_to;

    if(eventinfo->device) {
        selfimpot_from = Determine_self_coulomb(node1, eventinfo);
        selfimpot_to = Determine_self_coulomb(node2, eventinfo);   
        
        lr_coulomb_from = charge*Determine_lr_coulomb(node1, longrange, eventinfo);
        lr_coulomb_to = charge*Determine_lr_coulomb(node2, longrange, eventinfo);
        
        init_energy = static_node_energy_from + selfimpot_from + from_event_energy + sr_coulomb_from + lr_coulomb_from;
        final_energy = static_node_energy_to + selfimpot_to + to_event_energy + sr_coulomb_to + lr_coulomb_to;
    }
    else {
        init_energy = static_node_energy_from + from_event_energy + sr_coulomb_from;
        final_energy = static_node_energy_to + to_event_energy + sr_coulomb_to;
    }

    double energycontrib = 0.0;

    if (eventinfo->formalism == "Miller") {
        if(_final_type == Blocking) {
            _energyfactor = 0.0; // Keep this here for eventual simulation of bipolaron formation for example
        }
        else {
            energycontrib = final_energy - init_energy - charge*(eventinfo->efield_x*distancevector.x()+eventinfo->efield_y*distancevector.y()+eventinfo->efield_z*distancevector.z());
            if (energycontrib>0.0) {
                _energyfactor = exp(-1.0*energycontrib/(kB*eventinfo->temperature));
            }
            else {
                _energyfactor = 1.0;
            }
        }
    }
    else if (eventinfo->formalism == "Marcus") {
        if(_final_type == Blocking) {
            _energyfactor = 0.0; // Keep this here for eventual simulation of bipolaron formation for example
        }
        else {
            energycontrib = final_energy - init_energy - charge*(eventinfo->efield_x*distancevector.x()+eventinfo->efield_y*distancevector.y()+eventinfo->efield_z*distancevector.z());
            _energyfactor = exp(-1.0*(energycontrib+Reorg)*(energycontrib+Reorg)/(4*Reorg*kB*eventinfo->temperature));
        }       
    }
    _rate = prefactor*_transferfactor*_energyfactor;
    //if(node1->type() == (int) LeftElectrodeNode) std::cout << _rate << " " << _energyfactor << " " << final_energy << " " << static_node_energy_to << " " << _transferfactor << " " << distancevector << " " << (2*Pi/hbar)*(1.0/sqrt(4*Pi*Reorg*kB*eventinfo->temperature)) << " " << node1->position() << " " << node2->position() << " " << endl;

    //if(node1->type() == (int) RightElectrodeNode && _rate > 1.0) std::cout << _rate << " " << _energyfactor << " " << final_energy << " " << static_node_energy_to << " " << _transferfactor << " " << distancevector << " " << (2*Pi/hbar)*(1.0/sqrt(4*Pi*Reorg*kB*eventinfo->temperature)) << " " << node1->position() << " " << node2->position() << " " << endl;

}

void Event::Set_event(Link* link,int carrier_type, StateReservoir* state, Longrange* longrange, Eventinfo* eventinfo) {
   
    _link = link;
    Node* node1 = link->node1();
    Node* node2 = link->node2();
    
    _carrier_type = carrier_type;
    _init_type = Determine_init_event_type(node1);
    if (node2->occ() == -1) {_final_type = Determine_final_event_type(node1, node2);}
    else                    {_final_type = Determine_final_event_type(carrier_type, state->GetCarrier(node2->occ())->type(), node1, node2, eventinfo);}    

    _action_pair = (int) Other;
//    if(_action_pair != (int) Transfer) {
        _action_node1 = Determine_action_flag_node1();
        _action_node2 = Determine_action_flag_node2();
//    }
//    else {
//        _action_node1 = (int) None;
//        _action_node2 = (int) None;
//    }    
    
    Determine_rate(state, longrange, eventinfo);
}

}} 

#endif

