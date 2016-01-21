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

#ifndef __VOTCA_KMC_CARRIER_H_
#define __VOTCA_KMC_CARRIER_H_

#include <votca/xtp/node.h>

namespace votca { namespace xtp {

enum CarrierType{ Reservoir, Electron, Hole, Singlet, Triplet};    
    
class Carrier {
public:

    Carrier(int id){
        _carrier_id = id;
    };
    
    /// carrier_id
    const int &id() const { return _carrier_id; }
    /// carrier node
    Node* &node() { return _carrier_node; }    
    /// carrier type
    int &type() { return _carrier_type; } 
    /// travelled distance
    const votca::tools::vec &distance() const { return _carrier_distance; } 
    
    /// is carrier in box or not?
    bool &inbox() { return _in_sim_box; }
    
    /// set carrier node
    void SetCarrierNode(Node* carnode) { _carrier_node = carnode; }
    /// set travelled distance
    void SetDistance(votca::tools::vec distance) { _carrier_distance = distance; }
    void IncDistance(votca::tools::vec distance) { _carrier_distance += distance; }
    /// set "inbox/outbox" status
    void SetInBox(bool in_sim_box) {_in_sim_box = in_sim_box;}
    /// set carrier type
    void SetCarrierType(int carrier_type) {_carrier_type = carrier_type;}
    
    bool isElectron() { if(_carrier_type == Electron) {return true;} else {return false;}}
    bool isHole() { if(_carrier_type == Hole) {return true;} else {return false;}}
    bool isSinglet() { if(_carrier_type == Singlet) {return true;} else {return false;}}
    bool isTriplet() { if(_carrier_type == Triplet) {return true;} else {return false;}}
    
private:
    
    int _carrier_id;
    Node* _carrier_node;
    int _carrier_type;
    votca::tools::vec _carrier_distance;
    
    bool _in_sim_box;
};

}} 

#endif

