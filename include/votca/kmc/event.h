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
#include <votca/kmc/graphdevice.h>
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/state.h>
#include <votca/kmc/longrange.h>

namespace votca { namespace kmc {
  
using namespace std;

enum Final_Event_Type {TransferTo, Collection, CollectiontoNode, Recombination, Blocking, Notinbox};
enum Init_Event_Type {Injection, InjectionfromNode, TransferFrom, Notinboxfrom};
enum Action{Add, Remove, None };

class Event {
    
public:
    
    Event(int id, int type) {
        _id = id;
        _init_type = (int) Notinboxfrom;
        _final_type = type;
        _rate = 0.0;
    }

    Event(int id, Link* link, int carrier_type, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo){
        _id = id;
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
    
    void Set_event(Link* link, int carrier_type, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
    /// Determine rate
    void Determine_rate(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo);
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
    inline int Determine_final_event_type(int carrier_type1, int carrier_type2, Node* node1, Node* node2);
    /// Determine non injection event type in case only one carrier is on a link
    inline int Determine_final_event_type(Node* node1, Node* node2);
    /// Determine action flag for node 1
    inline int Determine_action_flag_node1();
    /// Determine action flag for node 2
    inline int Determine_action_flag_node2();

    /// Determine short range coulomb potential at node from which hop occur
    inline double Determine_from_sr_coulomb(Node* node1, StateDevice* state, Eventinfo* eventinfo); 
    
    /// Determine short range coulomb potential at node to which hop occur
    inline double Determine_to_sr_coulomb(Node* node1, StateDevice* state, Eventinfo* eventinfo);
    
    /// Determine long range coulomb potential for given node (equal to 0 in case of electrode node)
    inline double Determine_lr_coulomb(Node* node, Longrange* longrange);
    
protected:
                         
    Link* _link;
    int _carrier_type;
    
    double _rate;
    
    int _final_type;
    int _init_type;
    int _id;
    
    int _action_node1;
    int _action_node2;
    
    int _layer_node1;
    int _layer_node2;
    
    double _injection_potential;
    
};

inline int Event::Determine_final_event_type(Node* node1, Node* node2) {
    if (node2->type() == (int) NormalNode){                                                                  
        if(!dynamic_cast<NodeDevice*>(node2)->injectable()) return (int) TransferTo;
        if(dynamic_cast<NodeDevice*>(node2)->injectable()) return (int) CollectiontoNode;
    }
    else if (node2->type() == (int) LeftElectrodeNode || node2->type() == (int) RightElectrodeNode){         return (int) Collection;} // Collection at electrode
}

inline int Event::Determine_final_event_type(int carrier_type1, int carrier_type2, Node* node1, Node* node2) {
    if (((carrier_type1 == (int) Electron) && (carrier_type2 == (int) Electron)) || ((carrier_type1 == (int) Hole) && (carrier_type2 == (int) Hole))){          return (int) Blocking;     } // Blocking
    else if (((carrier_type1 == (int) Electron) && (carrier_type2 == (int) Hole))     || ((carrier_type1 == (int) Hole) && (carrier_type2 == (int) Electron))){ return (int) Recombination;} // Recombination
}

inline int Event::Determine_init_event_type(Node* node1) {
    if(node1->type() == (int) NormalNode){                                                                   
        if(!dynamic_cast<NodeDevice*>(node1)->injectable()) return (int) TransferFrom;
        if(dynamic_cast<NodeDevice*>(node1)->injectable()) return (int) InjectionfromNode;
    }
    else if((node1->type() == (int) LeftElectrodeNode) || (node1->type() == (int) RightElectrodeNode)){      return (int) Injection;   }
}

inline int Event::Determine_action_flag_node1() {
    int action_node1;
    if(_init_type == Injection)           {action_node1 = (int) None;   } // injection
    else if(_init_type == TransferFrom)   {action_node1 = (int) Remove; } // transfer
    return action_node1;    
}

inline int Event::Determine_action_flag_node2() {
    int action_node2;
    if(_final_type == TransferTo)         {action_node2 = (int) Add;    } // transfer
    else if(_final_type == Collection)    {action_node2 = (int) None;   } // collection
    else if(_final_type == Recombination) {action_node2 = (int) Remove; } // recombination
    return action_node2;    
}

inline double Event::Determine_from_sr_coulomb(Node* node, StateDevice* state, Eventinfo* eventinfo) {

    double coulomb_from;
    if(_init_type == Injection)   {coulomb_from = 0.0;                                                                           }
    else                          {coulomb_from = eventinfo->coulomb_strength*state->GetCarrier(node->occ())->on_site_coulomb(); }
    return coulomb_from;
}

inline double Event::Determine_to_sr_coulomb(Node* node, StateDevice* state, Eventinfo* eventinfo) {

    double coulomb_to;
    if(_final_type == Collection) {coulomb_to = 0.0;                                                                                       }
    else if(_init_type != Injection) {
        coulomb_to = eventinfo->coulomb_strength*state->GetCarrier(node->occ())->to_site_coulomb(_link->id());}
    else if(_init_type == Injection) {coulomb_to = eventinfo->coulomb_strength*_injection_potential;}
    return coulomb_to;
}

inline double Event::Determine_lr_coulomb(Node* node, Longrange* longrange) {
    double lr_coulomb;
    if(node->type() == (int) NormalNode) {lr_coulomb = longrange->Get_cached_longrange(dynamic_cast<NodeDevice*>(node)->layer());}
    else {lr_coulomb = 0.0;} //potential at electrodes = 0.0;
    return lr_coulomb;
}

void Event::Determine_rate(StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
    
    Node* node1 = _link->node1();
    Node* node2 = _link->node2();
    
    double prefactor = 1.0; // total prefactor
    double charge;
    double static_node_energy_from;
    double static_node_energy_to;

    const double hbar = 6.58211928E-16; // eV*s
    const double Pi   = 3.14159265358979323846264338327950288419716939937510;
    const double kB   = 8.617332478E-5; // ev/K     
    
    if(_carrier_type == (int) Electron) {
        charge = -1.0;
        prefactor = prefactor*(eventinfo->electron_prefactor);
        static_node_energy_from = dynamic_cast<NodeSQL*>(node1)->eCation() + dynamic_cast<NodeSQL*>(node1)->UcCnNe();
        static_node_energy_to = dynamic_cast<NodeSQL*>(node2)->eCation() + dynamic_cast<NodeSQL*>(node2)->UcCnNe();
    }
    else if(_carrier_type == (int) Hole) {
        charge = 1.0;
        prefactor = prefactor*(eventinfo->hole_prefactor);
        static_node_energy_from = dynamic_cast<NodeSQL*>(node1)->eAnion() + dynamic_cast<NodeSQL*>(node1)->UcCnNh();
        static_node_energy_to = dynamic_cast<NodeSQL*>(node2)->eAnion() + dynamic_cast<NodeSQL*>(node2)->UcCnNh();
        if(_init_type == Injection) static_node_energy_from += 21.5;
        if(_final_type == Collection) static_node_energy_to += 21.5;
    }
    
    //first transfer integrals
   
    double transferfactor = 1.0;
    votca::tools::vec distancevector = _link->r12();
    double Reorg;
    double Jeff2;
    
    if (eventinfo->formalism == "Miller") {
        double distance = abs(distancevector);
        transferfactor = exp(-1.0*eventinfo->alpha*distance);
    }
    else if((eventinfo->formalism == "Marcus")&&(_init_type==Injection||_final_type==Collection)) {
         double distance = abs(distancevector);
         transferfactor = exp(-2.0*eventinfo->alpha*distance);
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
        transferfactor = (2*Pi/hbar)*(Jeff2/sqrt(4*Pi*Reorg*kB*eventinfo->temperature));
    }
  
    //second, calculate boltzmann factor

    double init_energy;
    double final_energy;
    double from_event_energy = 0.0;
    double to_event_energy = 0.0;

    if(_init_type == Injection)           { from_event_energy -= eventinfo->injection_barrier; prefactor *= eventinfo->injection_prefactor;    } // injection
    else if(_init_type == TransferFrom)   {                                                                                                    } // transfer

    if(_final_type == TransferTo)         {                                                                                                    } // transfer
    else if(_final_type == Collection)    { to_event_energy   -= eventinfo->injection_barrier; prefactor *= eventinfo->collection_prefactor;   } // collection
    else if(_final_type == Recombination) { to_event_energy   -= eventinfo->binding_energy;    prefactor *= eventinfo->recombination_prefactor;} // recombination

    double sr_coulomb_from = Determine_from_sr_coulomb(node1, state, eventinfo);
    double sr_coulomb_to = Determine_to_sr_coulomb(node1, state, eventinfo);

    if(eventinfo->device) {
        double selfimpot_from = dynamic_cast<NodeDevice*>(node1)->self_image();
        double selfimpot_to = dynamic_cast<NodeDevice*>(node2)->self_image();    

        double lr_coulomb_from = eventinfo->coulomb_strength*charge*Determine_lr_coulomb(node1, longrange);
        double lr_coulomb_to = eventinfo->coulomb_strength*charge*Determine_lr_coulomb(node2, longrange);

        init_energy = static_node_energy_from + selfimpot_from + from_event_energy + sr_coulomb_from + lr_coulomb_from;
        final_energy = static_node_energy_to + selfimpot_to + to_event_energy + sr_coulomb_to + sr_coulomb_to;
    }
    else {
        init_energy = static_node_energy_from + from_event_energy + sr_coulomb_from;
        final_energy = static_node_energy_to + to_event_energy + sr_coulomb_to;        
    }

    double energycontrib;
    double energyfactor;

    if ((eventinfo->formalism == "Miller")||(_init_type==Injection||_final_type==Collection)) {
        if(_final_type == Blocking) {
            energyfactor = 0.0; // Keep this here for eventual simulation of bipolaron formation for example
        }
        else {
            energycontrib = final_energy - init_energy -charge*(eventinfo->efield_x*distancevector.x()+eventinfo->efield_y*distancevector.y()+eventinfo->efield_z*distancevector.z());
            if (energycontrib>0.0) {
                energyfactor = exp(-1.0*energycontrib/(kB*eventinfo->temperature));
            }
            else {
                energyfactor = 1.0;
            }
        }
    }
    else if (eventinfo->formalism == "Marcus") {
        if(_final_type == Blocking) {
            energyfactor = 0.0; // Keep this here for eventual simulation of bipolaron formation for example
        }
        else {
            energycontrib = final_energy - init_energy -charge*(eventinfo->efield_x*distancevector.x()+eventinfo->efield_y*distancevector.y()+eventinfo->efield_z*distancevector.z());
            energyfactor = exp(-1.0*(energycontrib+Reorg)*(energycontrib+Reorg)/(4*Reorg*kB*eventinfo->temperature));
        }       
    }

    _rate = prefactor*transferfactor*energyfactor;
//        std::cout << prefactor << " " << transferfactor << " " << energyfactor << " " << dynamic_cast<LinkSQL*>(_link)->Jeff2h() << " " << dynamic_cast<LinkSQL*>(_link)->lOh() << " " << Reorg << endl;
}

void Event::Set_event(Link* link,int carrier_type, StateDevice* state, Longrange* longrange, Eventinfo* eventinfo) {
    
    _link = link;
    Node* node1 = link->node1();
    Node* node2 = link->node2();
    
    _carrier_type = carrier_type;
    _init_type = Determine_init_event_type(node1);
    
    if (node2->occ() == -1) {_final_type = Determine_final_event_type(node1, node2);}
    else                    {_final_type = Determine_final_event_type(carrier_type, state->GetCarrier(node2->occ())->type(), node1, node2);}    

    _action_node1 = Determine_action_flag_node1();
    _action_node2 = Determine_action_flag_node2();
    
    Determine_rate( state, longrange, eventinfo);

}

}} 

#endif

