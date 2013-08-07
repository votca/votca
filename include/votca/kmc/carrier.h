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

#ifndef __VOTCA_KMC_CARRIER_H_
#define __VOTCA_KMC_CARRIER_H_

#include <votca/tools/vec.h>
#include <votca/kmc/node.h>

namespace votca { namespace kmc {
  
using namespace std;

//rewrite this as being a property of node (meaning, occupation and its actual 

class Carrier {
public:

    const int &getCarrierID() const { return _carrierID; }
    const CarrierType &getType() const { return _type; }
    const Node *getNode() const { return _node; }
    const vec &getDistance() const { return _distance; }
    
    void setType(CarrierType type) {
        _type = type;
    }
    
    void setNode(Node* position) {
        _node = position;
    }
    
    void setCarrierID(int Id) {
        _carrierID = Id
    }
    
    void 
    
private:
    CarrierType _type;
    vec _distance;
    Node* _node;
    int _carrierID;    
};

}} 

#endif

