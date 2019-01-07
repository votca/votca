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

bool compareChainWithChains_(vector<int> chain, vector<vector<int>> chains){
  bool match = false;
  for( auto chain2 : chains){
    if(chain2.size()==chain.size()){
      match = true;
      for(size_t index =0; index<chain.size();++index){
        if(chain2.at(index)!=chain.at(index)){
          match = false;
          break;
        }
      }// Cycle vertices in each chain
      if(match) return true;
    }// Chains same size
  }
  return false;
}


ReducedGraph::ReducedGraph(std::vector<ReducedEdge> reduced_edges){
  std::vector<Edge> edges;
  for(auto reduced_edge : reduced_edges){
    Edge ed(reduced_edge.getEndPoint1(),reduced_edge.getEndPoint2());
    auto temp_chain = reduced_edge.getChain();
    bool match = false;
    // Cycle the chains
    if(expanded_edges_.count(ed)){
      auto chains = expanded_edges_[ed];
      match = compareChainWithChains_(temp_chain, chains);
    }
     
    if(!match){
      edges.push_back(ed);
      expanded_edges_[ed].push_back(temp_chain);
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
  for(auto edge_and_chains : expanded_edges_ ){
    for(int edge_count=0; edge_count<edge_and_chains.second.size();++edge_count){
      edges.push_back(edge_and_chains.first);
    }
  }
  return edges;
}

ostream& operator<<(ostream& os, const ReducedGraph g){
  os << "Graph" << endl;
  for(auto p_gn : g.nodes_){
    os << "Node " << p_gn.first << endl;
    os << p_gn.second << endl;
  }

  os << endl;
  os << g.edge_container_ << endl;
  
  os << endl;
  os << "Expanded Edge Chains" << endl;
  
  for(auto edge_pr : g.expanded_edges_){
    for(auto chain : edge_pr.second){
      for(auto vertex : chain ){
        os << vertex << " ";
      }
      os << endl;
    }
    os << endl;
  }

  return os;
}



}
}
