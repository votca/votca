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
        for (it = _nodes.begin(); it != _nodes.end(); it++ ) delete *it;
    };   
    
    /// Add a node to the Graph
    void AddNode(Node* node) { _nodes.push_back(node); }
    
    void Print(std::ostream& out){
        std::vector<Node*>::iterator it;
          for (it = _nodes.begin(); it != _nodes.end(); it++ ) (*it)->Print(out);    
    }
    
    /// initialize nodes and links
    virtual void Initialize() = 0;
    
private:
    std::vector<Node*> _nodes;
    std::vector<Link*> _links;
    
};


}}



#endif

