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

#include <votca/tools/reducedgraph.h>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;

ReducedGraph::ReducedGraph(std::vector<ReducedEdge> reduced_edges){
  std::vector<Edge> edges;
  for(auto reduced_edge : reduced_edges){
    Edge ed(reduced_edge.getV1(),reduced_edge.getV2());
    auto temp_chain = reduced_edge.getChain();
    bool match = false;
    if(edge_count_.count(ed)){
      match = false;
      for(auto edge_chain_ptr : edge_count_[ed]){
        if((*edge_chain_ptr).second.size()==temp_chain.size()){
          for(size_t i=0; i<temp_chain.size();++i){
            if(temp_chain.at(i)!=(*edge_chain_ptr).second.at(i)){
              break;
            }
            match=true;
          }
        }
        if(match) break;        
      }
    }
    
    if(!match){
      edges.push_back(ed);
      expanded_edges_.push_back(pair<Edge,vector<int>>(ed, temp_chain));
      edge_count_[ed].push_back(&expanded_edges_.back());
    }
  }
  edge_container_ = EdgeContainer(edges);
  auto vertices = edge_container_.getVertices();
  for(auto vertex : vertices ){
    GraphNode gn;
    nodes_[vertex] = gn;  
  }
}

ReducedGraph::ReducedGraph(const ReducedGraph& g) {
  this->edge_container_ = g.edge_container_;
  this->expanded_edges_ = g.expanded_edges_;
  for(auto pr : g.nodes_ ){
    this->nodes_[pr.first] = pr.second;
  }
  this->id_ = g.id_;
}

ReducedGraph& ReducedGraph::operator=(const ReducedGraph& g) {
  this->expanded_edges_ = g.expanded_edges_;
  this->edge_container_ = g.edge_container_;
  for( auto pr : g.nodes_ ){
    this->nodes_[pr.first] = pr.second;
  }
  this->id_ = g.id_;
  return *this;
}

ReducedGraph& ReducedGraph::operator=(ReducedGraph&& g) {
  this->expanded_edges_ = move(g.expanded_edges_);
  this->edge_container_ = move(g.edge_container_);
  this->nodes_ = move(g.nodes_);
  this->id_ = move(g.id_);
  return *this;
}

vector<Edge> ReducedGraph::getEdges(){
  vector<Edge> edges;
  for(auto edge_and_vec : edge_count_ ){
    for(int edge_count=0; edge_count<edge_and_vec.second.size();++edge_count){
      edges.push_back(edge_and_vec.first);
    }
  }
  return edges;
}
}
}
