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

namespace votca { namespace kmc {
  
using namespace std;

enum Final_Event_Type {TransferTo, Collection, Recombination, Blocking, Notinbox};
enum Init_Event_Type {Injection, TransferFrom};
enum Action{Add, Remove, None };

class Event {
    
public:

    Event(Link* link, Carrier* carrier, Eventinfo* eventinfo, StateDevice<GraphDevice<GraphSQL, NodeSQL, LinkSQL> >* state){
        _link = link;
        Set_event(link, carrier, eventinfo, state);
    };
 
    void Set_event(Link* link, Carrier* carrier, Eventinfo* eventinfo, StateDevice<GraphDevice<GraphSQL, NodeSQL, LinkSQL> >* state);

    /// Determine initial event type
    int Determine_init_event_type(Node* node1);
    /// Determine non injection event type in case two carriers are on linking nodes
    int Determine_final_event_type(Carrier* carrier1, Carrier* carrier2, Node* node1, Node* node2);
    /// Determine non injection event type in case only one carrier is on a link
    int Determine_final_event_type(Node* node1, Node* node2);

    
 //   void Set_injection_event(DNode* electrode, int injectnode_ID, CarrierType carrier_type,
 //                                double from_longrange, double to_longrange, Globaleventinfo* globevent);
 //   void Set_non_injection_event(vector<DNode*> nodes, Carrier* carrier, int jump_ID,
 //                                double from_longrange, double to_longrange, Globaleventinfo* globevent);
    
protected:
 //   From_step_event Determine_non_injection_from_event_type(Carrier* carrier);
 //   To_step_event Determine_non_injection_to_event_type(Carrier* carrier, int jumpID, DNode* carriernode);
 //   To_step_event Determine_injection_to_event_type(CarrierType carrier_type, DNode* electrode, int inject_nodeID);

 //   double Compute_event_rate(DNode* fromnode, int jump_ID, CarrierType carrier_type,
 //                           From_step_event from_event_type, To_step_event to_event_type,
//                            double from_shortrange, double to_shortrange, double from_longrange, double to_longrange,
//                            Globaleventinfo* globaleventinfo);
    
    Link* _link;
    Carrier* _carrier;
    
    double _rate;
    int _final_type;
    int _init_type;
    
    int _action_node1;
    int _action_node2;
    
};

int Event::Determine_final_event_type(Node* node1, Node* node2) {
    if (node2->type() == (int) NormalNode)                                                                  return (int) TransferTo; // Transfer to empty node
    if (node2->type() == (int) LeftElectrodeNode || node2->type() == (int) RightElectrodeNode)              return (int) Collection; // Collection at electrode
}

int Event::Determine_final_event_type(Carrier* carrier1, Carrier* carrier2, Node* node1, Node* node2) {
    if ((carrier1->isElectron() && carrier2->isElectron())||(carrier1->isHole() && carrier2->isHole()))     return (int) Blocking; // Blocking
    if ((carrier1->isElectron() && carrier2->isHole())    ||(carrier1->isHole() && carrier2->isElectron())) return (int) Recombination; // Recombination
}

int Event::Determine_init_event_type(Node* node1) {
    if(node1->type() == (int) NormalNode)                                                                   return (int) TransferFrom;
    if((node1->type() == (int) LeftElectrodeNode) || (node1->type() == (int) RightElectrodeNode))           return (int) Injection;
}

void Event::Set_event(Link* link, Carrier* carrier, Eventinfo* eventinfo, StateDevice<GraphDevice<GraphSQL, NodeSQL, LinkSQL> >* state) {
    
    if(!carrier->inbox()){ _final_type = (int) Notinbox; _rate = 0.0; }
    else {
    
        Node* node1 = link->node1();
        Node* node2 = link->node2();

        _init_type = Determine_init_event_type(node1);
        
        if (node2->occ() == -1) _final_type = Determine_final_event_type(node1, node2);
        else                    _final_type = Determine_final_event_type(carrier, state->GetCarrier(node2->occ()), node1, node2);    

        double prefactor = 1.0; // total prefactor

        double charge;
        double static_node_energy_from;
        double static_node_energy_to;

        if(carrier->isElectron()) {
            charge = 1.0;
            prefactor *= eventinfo->electron_prefactor;
            static_node_energy_from = dynamic_cast<NodeSQL*>(node1)->eCation() + dynamic_cast<NodeSQL*>(node1)->ucCnNe();
            static_node_energy_to = dynamic_cast<NodeSQL*>(node2)->eCation() + dynamic_cast<NodeSQL*>(node2)->ucCnNe();
        }
        else if(carrier->isHole()) {
            charge = -1.0;
            prefactor *= eventinfo->hole_prefactor;
            static_node_energy_from = dynamic_cast<NodeSQL*>(node1)->eAnion() + dynamic_cast<NodeSQL*>(node1)->ucCnNh();
            static_node_energy_to = dynamic_cast<NodeSQL*>(node2)->eAnion() + dynamic_cast<NodeSQL*>(node2)->ucCnNh();
        }

        
        //first calculate quantum mechanical wavefunction overlap
        votca::tools::vec distancevector = link->r12();
        double distance = abs(distancevector);

        double distancefactor = exp(-1.0*eventinfo->alpha*distance);

        //second, calculate boltzmann factor (Coulomb interaction still have to be implemented)

        double init_energy;
        double final_energy;
    //    double selfimpot_from = fromnode->self_image_potential;
    //    double selfimpot_to = jumptonode->self_image_potential;
        double from_event_energy = 0.0;
        double to_event_energy = 0.0;

    //    if(from_event_type == Injection) {
    //        fromtype_energy -= globevent->injection_barrier;
    //        prefactor *= globevent->injection_prefactor;
    //    }
        if(_init_type == Injection)           { from_event_energy -= eventinfo->injection_barrier; prefactor *= eventinfo->injection_prefactor;     _action_node1 = (int) None;   } // injection
        else if(_init_type == TransferFrom)   {                                                                                                     _action_node1 = (int) Remove; } // transfer
        
        if(_final_type == TransferTo)         {                                                                                                     _action_node2 = (int) Add;    } // transfer
        else if(_final_type == Collection)    { to_event_energy   -= eventinfo->injection_barrier; prefactor *= eventinfo->collection_prefactor;    _action_node2 = (int) None;   } // collection
        else if(_final_type == Recombination) { to_event_energy   -= eventinfo->binding_energy;    prefactor *= eventinfo->recombination_prefactor; _action_node2 = (int) Remove; } // recombination

    //    double coulomb_from;
    //    double coulomb_to;

    //    if ((from_event_type != Injection)) { 
    //        coulomb_from = globevent->coulomb_strength*(from_shortrange+charge*from_longrange);
    //        coulomb_to = globevent->coulomb_strength*(to_shortrange+charge*to_longrange);
    //    }
    //    else if(from_event_type == Injection){
    //        coulomb_from = 0.0;
    //        coulomb_to = charge*globevent->coulomb_strength*(jumptonode->injection_potential+to_longrange);
    //    }

    //    init_energy = static_node_energy_from + selfimpot_from + coulomb_from;
        init_energy = static_node_energy_from;
    //    final_energy = static_node_energy_to + selfimpot_to + coulomb_to;
        final_energy = static_node_energy_to;

        double energycontrib;
        double energyfactor;

        if (eventinfo->formalism == "Miller") {
            if(_final_type == Blocking) {
                energyfactor = 0.0; // Keep this here for eventual simulation of bipolaron formation for example
            }
            else {
                energycontrib = final_energy - init_energy -charge*(eventinfo->efield_x*distancevector.x()+eventinfo->efield_y*distancevector.y()+eventinfo->efield_z*distancevector.z());
                if (energycontrib>0.0) {
                    energyfactor = exp(-1.0*eventinfo->beta*energycontrib);
                }
                else {
                    energyfactor = 1.0;
                }
            }
        }

        _rate = prefactor*distancefactor*energyfactor;
   }
}

/*void Event::Set_injection_event(DNode* electrode, int injectnode_ID, CarrierType carrier_type,
                              double from_longrange, double to_longrange, Globaleventinfo* globevent) {
    
    fromtype = Injection;
    totype = Determine_injection_to_event_type(carrier_type, electrode, injectnode_ID);
    rate = Compute_event_rate(electrode, injectnode_ID, carrier_type, fromtype, totype,
                              0, 0.0, from_longrange, to_longrange, globevent);
       
}*/
/*

void Event::Set_non_injection_event(vector<DNode*> nodes, Carrier* carrier, int jump_ID,
                                 double from_longrange, double to_longrange, Globaleventinfo* globaleventinfo) {
    
    fromtype = Determine_non_injection_from_event_type(carrier);
    totype = Determine_non_injection_to_event_type(carrier, jump_ID, nodes[carrier->carrier_node_ID]);
    rate = Compute_event_rate(nodes[carrier->carrier_node_ID], jump_ID, carrier->carrier_type, fromtype, totype,
                                     carrier->srfrom, carrier->srto[jump_ID], from_longrange, to_longrange,
                                     globaleventinfo);    
}
*/



}} 

#endif

