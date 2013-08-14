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

typedef votca::tools::vec myvec;

namespace votca { namespace kmc {
  
using namespace std;

enum CarrierType{ Electron, Hole, Exciton, Void };

class Carrier {
public:

    const int &getCarrierID() const { return _carrierID; }
    const CarrierType &getCarrierType() const { return _carriertype; }
    const Node* getNode() const { return _node; }
    const myvec &getDistance() const { return _distance; }
    
//    void setType(CarrierType type) {
//        _type = type;
//    }
    
//    void setCarrierID(int Id) {
//        _carrierID = Id;
//    }

    CarrierType _carriertype;
    myvec _distance;
    Node* _node;
    int _carrierID;    
};

}} 

#endif

