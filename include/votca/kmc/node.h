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

#ifndef __VOTCA_KMC_NODE_H_
#define	__VOTCA_KMC_NODE_H_

#include <vector>
#include <votca/tools/vec.h>
#include "carrier.h"

typedef votca::tools::vec myvec;

namespace votca { namespace kmc {
  
using namespace std;

enum NodeType {Normal, LeftElectrode, RightElectrode};


class Node {
    
  public:
        
    double &x() { return _position.x(); }   
    double &y() { return _position.y(); }
    double &z() { return _position.z(); }
        
    void setElPair(Node* finalnode, myvec distance);
    void setHoPair(Node* finalnode, myvec distance);
    void setExcitonPair(Node* finalnode, myvec distance);
        
    void setElTransfers(double rate12e, double Jeff2e, double reorg_oute);
    void setHoTransfers(double rate12h, double Jeff2h, double reorg_outh);
        
    int _node_id;
    myvec _position;
        
    // Parameters needed for Marcus Formalism 
    double _reorg_intorig_el;
    double _reorg_intorig_ho;
    double _reorg_intdest_el; 
    double _reorg_intdest_ho;
    double _eAnion;
    double _eCation;
    double _eNeutral;
    double _internalenergy_el;
    double _internalenergy_ho;        

    // General energy parameters
    double _static_electron_energy;
    double _static_hole_energy;        

    // Define neighbouring nodes
    vector<Node*> _el_pair_nodes;
    vector<Node*> _ho_pair_nodes;
    vector<Node*> _exciton_pair_nodes;
    vector<myvec> _el_pair_distances;  
    vector<myvec> _ho_pair_distances;
    vector<myvec> _exciton_pair_distances;
    int _number_of_ho_pairs;
    int _number_of_el_pairs;
    int _number_of_exciton_pairs;
        
    // Define relevant rate information
    vector<double> _rate12_el;
    vector<double> _rate12_ho;
    vector<double> _Jeff2_el;
    vector<double> _Jeff2_ho;
    vector<double> _reorg_out_el;
    vector<double> _reorg_out_ho;

    // Define relevant device setup information
    bool _left_injectable;
    bool _right_injectable;
    int _left_injectorindex;
    int _right_injectorindex;
    int _injectorindex;
    NodeType _node_type;       

    int _hole_number;
    int _electron_number;
    Carrier* _carrier;
    
    int _layer_index;
        
    double _el_injection_potential;
    double _ho_injection_potential;
};

       
void Node::setElPair(Node* finalnode, myvec distance){
  _el_pair_nodes.push_back(finalnode);
  _el_pair_distances.push_back(distance);
}
        
void Node::setHoPair(Node* finalnode, myvec distance){
  _ho_pair_nodes.push_back(finalnode);
  _ho_pair_distances.push_back(distance);
}

void Node::setExcitonPair(Node* finalnode, myvec distance){
  _exciton_pair_nodes.push_back(finalnode);
  _exciton_pair_distances.push_back(distance);
}
        
void Node::setElTransfers(double rate12e, double Jeff2e, double reorg_oute) {
  _rate12_el.push_back(rate12e) ;
  _Jeff2_el.push_back(Jeff2e);
  _reorg_out_el.push_back(reorg_oute);
}

void Node::setHoTransfers(double rate12h, double Jeff2h, double reorg_outh) {
  _rate12_ho.push_back(rate12h) ;
  _Jeff2_ho.push_back(Jeff2h);
  _reorg_out_ho.push_back(reorg_outh);            
}

/*        void Node::setExcitonTransfers(double rate12exc, double Jeff2exc, double reorg_outexc) {
            
        }*/        
        
}} 

#endif

