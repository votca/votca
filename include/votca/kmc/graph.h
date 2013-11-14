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
#include <votca/kmc/nodesql.h>
#include <votca/kmc/linksql.h>

namespace votca { namespace kmc {

template<class TNode, class TLink>
class Graph {

public:
     Graph() {}
     
    ~Graph() {
        typename std::vector<TNode*>::iterator it;
        for (it = _nodes.begin(); it != _nodes.end(); it++ ) delete *it;
    } 
    
    /// Add a node to the Graph
    TNode* AddNode(int id, votca::tools::vec &position) { 
        TNode* node = new TNode(id, position);
        _nodes.push_back(node); 
        return node;
    } 

    void AddNode( TNode* node) { _nodes.push_back(node); }
    
    TNode* GetNode(int node_ID) {return _nodes[node_ID];}

    /// Add a node to the Graph
    TLink* AddLink( int id, TNode* node1, TNode* node2, votca::tools::vec r12) { 
        TLink* link = new TLink(id, node1, node2, r12);
        _links.push_back(link); 
        return link;
    }
    
    void AddLink( TLink* link) { _links.push_back(link); }
    
    void PrintNodes(std::ostream& out){
        typename std::vector<TNode*>::iterator it;
//          for (it = _nodes.begin(); it != _nodes.end(); it++ ) (*it)->Print(out);    
        for (it = _nodes.begin(); it != _nodes.end(); it++ ) out << (*it)->id() << " " << (*it)->position() << " " <<
               (*it)->UnCnNe() << " " << (*it)->UnCnNh() << " " << (*it)->UcNcCe() << " " << (*it)->UcNcCh() << " " <<
               (*it)->eAnion() << " " << (*it)->eNeutral() << " " << (*it)->eCation() << " " <<
               (*it)->ucCnNe() << " " << (*it)->ucCnNh() << endl;    
    }

    void PrintLinks(std::ostream& out){
        typename std::vector<TLink*>::iterator it;
//          for (it = _nodes.begin(); it != _nodes.end(); it++ ) (*it)->Print(out);    
        for (it = _links.begin(); it != _links.end(); it++ ) out << (*it)->id() << " " << 
               (*it)->r12() << " " << (*it)->rate12e() << " " << (*it)->rate12h() << " " << (*it)->rate21e() << " " << (*it)->rate21h() << " " <<
               (*it)->Jeff2e() << " " << (*it)->Jeff2h() << " " << (*it)->lOe() << " " << (*it)->lOh() << endl;    
    }
    
    /// initialize nodes and links
    virtual void Initialize(string filename){;};

   
    
protected:

    std::vector<TNode*> _nodes;
    std::vector<TLink*> _links;     
};



}}



#endif

