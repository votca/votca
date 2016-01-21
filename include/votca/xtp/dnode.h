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

#ifndef __VOTCA_KMC_DNODE_H_
#define	__VOTCA_KMC_DNODE_H_

#include <vector>
#include <votca/tools/vec.h>
#include <votca/xtp/carrier.h>
#include <votca/xtp/node.h>

typedef votca::tools::vec myvec;

namespace votca { namespace xtp {
  
using namespace std;

enum Node_Type{Normal, LeftElectrode, RightElectrode};
/*
class DNode: public Node 
{
   
public:

    void setPair(DNode* pairing_node) {pairing_nodes.push_back(pairing_node);}
    void setStaticeventinfo(DNode* pairnode, myvec dr, double rate12e, double rate12h, double Jeff2e, double Jeff2h, double reorg_oute, double reorg_outh);    

    void removePair(int pairing_node_index);
    
    struct Static_event_info {
        DNode* pairnode;
        myvec distance; //distance vector from start to destination node
        double rate12e;
        double rate12h;
        double Jeff2e;
        double Jeff2h;
        double reorg_oute;
        double reorg_outh;
    };     

    int node_ID;
    Node_Type node_type;
    myvec node_position;
    vector<DNode*> pairing_nodes;
    vector<Static_event_info> static_event_info;
    vector<Carrier*> carriers_on_node;
    
    int layer_index;
    
    //static energies
    double reorg_intorig_hole;
    double reorg_intorig_electron;
    double reorg_intdest_hole;
    double reorg_intdest_electron;
        
    double eAnion;
    double eNeutral;
    double eCation;
        
    double internal_energy_electron;
    double internal_energy_hole;
        
    double static_electron_node_energy;
    double static_hole_node_energy;
    
    double self_image_potential;
    
    //for injection
    
    int left_injector_ID;
    int right_injector_ID;
    double injection_potential;
    
};

inline void DNode::setStaticeventinfo(DNode* pairnode, myvec dr, double rate12e, double rate12h, double Jeff2e, double Jeff2h, double reorg_oute, double reorg_outh) {
    Static_event_info newStatic;
    newStatic.pairnode = pairnode;
    newStatic.distance = dr;
    newStatic.rate12e = rate12e;
    newStatic.rate12h = rate12h;
    newStatic.Jeff2e = Jeff2e;
    newStatic.Jeff2h = Jeff2h;
    newStatic.reorg_oute = reorg_oute;
    newStatic.reorg_outh = reorg_outh;
    static_event_info.push_back(newStatic);
}

inline void DNode::removePair(int pairing_node_index) {
    pairing_nodes.erase(pairing_nodes.begin()+pairing_node_index);
    static_event_info.erase(static_event_info.begin()+pairing_node_index);
}
   */     
}} 

#endif

