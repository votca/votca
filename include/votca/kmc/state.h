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

#ifndef __VOTCA_KMC_STATE_H_
#define __VOTCA_KMC_STATE_H_

#include <votca/kmc/event.h>
#include <votca/kmc/carrier.h>
#include <votca/kmc/node.h>
//#include <votca/kmc/graph.h>
#include <votca/kmc/store.h>

// Text archive that defines boost::archive::text_oarchive
// and boost::archive::text_iarchive
//#include <boost/archive/text_iarchive.hpp>
//nclude <boost/archive/text_oarchive.hpp>

//#include <boost/serialization/list.hpp>

namespace votca { namespace kmc {
  
using namespace std;


class State {
public:
    
// container methods    
    
    void Grow(int growsize);
    int Buy();
    void Sell(int itemnr); //check whether this can work with direct pointer
    Carrier* Get_item(int itemnr);
    
//    void Clear(int totalnumberofnodes);
    void Save(){};
    void Load(){};

//   template<typename Archive> 
//   void serialize(Archive& ar, const unsigned version) {
//       //ar & _carriers; 
//       ;
//   }     
    
    Store<Carrier> _carriers;  
    
};

void State::Grow(int growsize){
    _carriers.Grow(growsize);
}

int State::Buy(){
    _carriers.Buy();
}

void State::Sell(int itemnr){
    _carriers.Sell(itemnr);
}

Carrier* State::Get_item(int itemnr){
    _carriers.Get_item(itemnr);
}
    
/*void State::Clear(int totalnumberofnodes) {
    
    for (int node_index=0;node_index<totalnumberofnodes;node_index++) {
        nodes[node_index]->setOccupation(Electron,0.0);
        nodes[node_index]->setOccupation(Hole,0.0);
        nodes[node_index]->setOccupation(Exciton,0.0);
        nodes[node_index]->setCarrierID(-1);

    }
}*/

}} 

#endif

