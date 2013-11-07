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

#ifndef __VOTCA_KMC_GRAPH_H_
#define __VOTCA_KMC_GRAPH_H_

#include <vector>
#include <votca/kmc/node.h>

namespace votca { namespace kmc {
  
class Graph {

public:
     Graph() {};
     
    ~Graph() {
        std::vector<Node*>::iterator it;
        for (it = nodes.begin(); it != nodes.end(); it++ ) delete *it;
    };   
    
    /// Add a node to the Graph
    void AddNode(Node* node) { nodes.push_back(node); }
    
    void Print(std::ostream& outstream){
        std::vector<Node*>::iterator it;
        /*for (it = nodes.begin(); it != nodes.begin()+10; it++ ) std::cout << (*it)->node_ID << " " << (*it)->node_position << " " 
                << (*it)->reorg_intorig_hole << " " << (*it) ->reorg_intorig_electron << " "
                << (*it)->reorg_intdest_hole << " " << (*it) ->reorg_intdest_electron << " "
                << (*it)->eAnion << " " << (*it) ->eNeutral << " " << (*it) ->eCation << " "
                << (*it)->internal_energy_electron << " " << (*it) ->internal_energy_hole << " "
                << (*it)->static_electron_node_energy << " " << (*it)->static_hole_node_energy                   
                << endl;   */
          for (it = nodes.begin(); it != nodes.end(); it++ ) std::cout << (*it)->node_ID << " " << (*it)->node_position << endl;    
    }
    
private:
    std::vector<Node*> nodes;
    
};


}}



#endif

