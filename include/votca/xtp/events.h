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

#ifndef __VOTCA_KMC_EVENTS_H_
#define __VOTCA_KMC_EVENTS_H_


#include <votca/xtp/graphkmc.h>
#include <votca/xtp/statereservoir.h>
#include <votca/xtp/event.h>
#include <votca/xtp/bsumtree.h>
#include <votca/xtp/longrange.h>
#include <votca/xtp/eventinfo.h>

namespace votca { namespace xtp {
  
using namespace std;

class Events {

public: 
    Events() {
    }
     
    ~Events() {
        std::vector<Event*>::iterator it;
        for (it = _non_injection_events.begin(); it != _non_injection_events.end(); it++ ) delete *it;
        for (it = _injection_events.begin(); it != _injection_events.end(); it++ ) delete *it;
    } 
    
    /// On execute method in device
    void On_execute(Event* event, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo);
    /// On execute method node wise
    void On_execute_node(Node* node, int action, int carrier_type, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo);
    /// On execute method pair wise
    void On_execute_pair(Link* link, int carrier_type, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo);
    
    /// Execute method in case carrier is added on node
    void Add_carrier(Node* node, int carrier_type, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo);

    /// Execute method in case carrier is removed from node
    void Remove_carrier(Node* node, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo);

    /// Recalculate rates of all events
    void Recompute_all_events(StateReservoir* state, Longrange* longrange,Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates,  Eventinfo* eventinfo);
    /// Recalculate rates of all non-injection events
    void Recompute_all_non_injection_events(StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Eventinfo* eventinfo);
    /// Recalculate rates of all injection events
    void Recompute_all_injection_events(StateReservoir* state, Longrange* longrange, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo);    

    /// Initialize event vectors
    void Initialize_device_eventvector(GraphKMC* graph, StateReservoir* state, Longrange* longrange, Eventinfo* eventinfo);
    
    void Initialize_bulk_eventvector(GraphKMC* graph, StateReservoir* state, Eventinfo* eventinfo);

    /// Initialize injection event vector
    void Initialize_injection_eventvector(int Event_counter, Node* electrode, int carrier_type, StateReservoir* state, Longrange* longrange, Eventinfo* eventinfo);
    /// Grow (and initialize) non-injection event vector
    void Grow_non_injection_eventvector(StateReservoir* state, Eventinfo* eventinfo);
    
    /// Initialize rates (after initialization of events)
    void Initialize_rates(Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo);

    /// Initialize mesh/potential and rates after placement of charges
    void Initialize_after_charge_placement(GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo);
    
    /// Initialize mesh for non-injection events
    void Init_non_injection_meshes(Eventinfo* eventinfo);
    /// Initialize mesh for injection events
    void Init_injection_meshes(StateReservoir* state, Eventinfo* eventinfo);
    /// Resize mesh
    vector< vector< vector <list<int> > > > Resize_mesh(int meshnr_x, int meshnr_y, int meshnr_z);
    /// Add id (of node or carrier) to mesh
    void Add_to_mesh(int ID, votca::tools::vec position, Eventinfo* eventinfo);
    /// Remove id (of node or carrier) from mesh
    void Remove_from_mesh(int ID,votca::tools::vec position, Eventinfo* eventinfo);

    /// Effect of adding/removing carrier to/from box on the coulomb potential and event rates
    void Effect_potential_and_rates(int action, CarrierBulk* carrier, Node* node, GraphKMC* graph, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Bsumtree* left_injection_rates, Bsumtree* right_injection_rates, Eventinfo* eventinfo);    

    /// Effect of adding/removing carrier to/from box on the coulomb potential and event rates (non-injection)
    void Effect_potential_and_non_injection_rates(int action, CarrierBulk* carrier1, Node* node, StateReservoir* state, Longrange* longrange, Bsumtree* non_injection_rates, Eventinfo* eventinfo);
    
    void Effect_non_injection_rates(int action, CarrierBulk* carrier, Node* node, StateReservoir* state, Bsumtree* non_injection_rates, Eventinfo* eventinfo);
    /// Effect of adding/removing carrier to/from box on the coulomb potential and event rates (injection)
    void Effect_injection_rates(int action, CarrierBulk* carrier, Node* node, Node* electrode, double dist_to_electrode, StateReservoir* state, Longrange* longrange, Bsumtree* injection_rates,  Eventinfo* eventinfo);
    
    /// Calculate shortrange coulomb potential
    /// startz is z-coordinate of carrier of which we want to calculate the coulomb potential
    /// dif is startz - z coordinate of node on which we want to know the coulomb potential
    double Compute_Coulomb_potential(double startz, votca::tools::vec dif, bool direct, votca::tools::vec sim_box_size, Eventinfo* eventinfo);

    /// Obtain non-injection event with id eventID
    Event* get_non_injection_event(int eventID) {return _non_injection_events[eventID];}
    /// Obtain injection event with id eventID
    Event* get_injection_event(int electrodeID, int eventID) {
        if(electrodeID == 0) {
            return _injection_events[eventID];
        }
        else if(electrodeID == 1) {
            return _injection_events[eventID + _total_left_injection_events];
        }
        return NULL;
    }
    
    double av_rate(int electrode);
    double av_inject_transferfactor(int electrode);
    double av_inject_energyfactor(int electrode);
    
private:

    vector<Event*> _non_injection_events;
    vector<Event*> _injection_events;
 
    int _total_non_injection_events;
    int _total_left_injection_events;
    int _total_right_injection_events;    
    
    double _meshsize_x; double _meshsize_y; double _meshsize_z; int _inject_meshnr_z;
    vector< vector< vector <list<int> > > > _non_injection_events_mesh;
    vector< vector< vector <list<int> > > > _left_injection_events_mesh;
    vector< vector< vector <list<int> > > > _right_injection_events_mesh;
};


}} 

#endif
