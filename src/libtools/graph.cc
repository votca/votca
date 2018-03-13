/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <string>
#include <votca/tools/graph.h>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;

void Graph::updateStructIds_(Graph& g){
  if(!this->structure_id_set_){
    this->calcStructureId_();
    this->structure_id_set_=true;
  }
  if(!g.structure_id_set_){
    g.calcStructureId_();
    g.structure_id_set_=true;
  }
}

bool Graph::operator!=(Graph& g) {
  updateStructIds_(g);
  return structure_id_.compare(g.structure_id_);
}

bool Graph::operator==( Graph& g) {
  return !(*(this)!=g);  
}

vector<pair<int, GraphNode>> Graph::getIsolatedNodes(void){
  vector<pair<int,GraphNode>> iso_nodes;
  for(auto node : nodes_ ){
    if(adj_list_.count(node.first)){
      if(adj_list_[node.first].size()==0){
        pair<int,GraphNode> pr(node.first,node.second);
        iso_nodes.push_back(pr);
      }
    }else{
      pair<int,GraphNode> pr(node.first,node.second);
      iso_nodes.push_back(pr);
    }
  }
  return iso_nodes;
}

vector<int> Graph::getVerticesMissingNodes(void){
  vector<int> missing;
  for(auto pr_v : adj_list_){
    if(nodes_.count(pr_v.first)==0){
      missing.push_back(pr_v.first);
    }
  }
  return missing;
}

GraphNode& Graph::Node(int vert){
  return nodes_[vert];
}

GraphNode Graph::getNode(int vert){
  return nodes_[vert];
}

vector<pair<int,GraphNode>> Graph::getNodes(void){
  vector<pair<int,GraphNode>> vec_nodes;
  for(auto pr_node : nodes_ ){
    vec_nodes.push_back(pr_node);
  } 
  return vec_nodes;
}

void Graph::calcStructureId_(){
  throw runtime_error("calcStructureId_ is not yet implemented");
  return;
}

ostream& operator<<(ostream& os, const Graph g){
  os << "Graph" << endl;
  for( auto p_gn : g.nodes_ ) {
    os << "Node " << p_gn.first << endl;
    os << p_gn.second << endl;
  }
  return os;
}
}
}

