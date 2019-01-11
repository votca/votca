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

#include <cassert>
#include <votca/tools/edge.h>
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

set<int> getAllVertices_(std::vector<ReducedEdge> reduced_edges){
  set<int> vertices;
  for( auto reduced_edge : reduced_edges){
    vector<int> chain = reduced_edge.getChain();
    for( auto vertex : chain){
      vertices.insert(vertex);
    }
  }
  return vertices;
}

ReducedGraph::ReducedGraph(std::vector<ReducedEdge> reduced_edges){

  auto vertices = getAllVertices_(reduced_edges);
  unordered_map<int,GraphNode> nodes;
  for( auto vertex : vertices  ){
    GraphNode gn;
    nodes[vertex] = gn;  
  }
  init_(reduced_edges,nodes);
}

ReducedGraph::ReducedGraph(std::vector<ReducedEdge> reduced_edges, unordered_map<int,GraphNode> nodes){
  auto vertices = getAllVertices_(reduced_edges);
  if(nodes.size()<vertices.size()){
    throw invalid_argument("The number of nodes passed into a reduced graph "
        "must be greater or equivalent to the number of vertices");
  }
  for(auto vertex : vertices){
    if(nodes.count(vertex)==0){
      throw invalid_argument("A vertex is missing its corresponding node.");
    }
  }
  init_(reduced_edges,nodes);
}

void ReducedGraph::init_(std::vector<ReducedEdge> reduced_edges, unordered_map<int,GraphNode> nodes){
  vector<Edge> edges;
  nodes_ = nodes;
  for(auto reduced_edge : reduced_edges){
    Edge ed(reduced_edge.getEndPoint1(),reduced_edge.getEndPoint2());
    edges.push_back(ed);
    if(expanded_edges_.count(ed)){
      bool match = compareChainWithChains_(reduced_edge.getChain(),expanded_edges_[ed]);
      if(!match){
        expanded_edges_[ed].push_back(reduced_edge.getChain());
      }
    }else{
      expanded_edges_[ed].push_back(reduced_edge.getChain());
    }
  }
  edge_container_ = EdgeContainer(edges);
  
  auto edges_test = edge_container_.getEdges();
  cout << "Size of edges " << edges_test.size() << endl;
  calcId_();
}

ReducedGraph::ReducedGraph(const ReducedGraph& g) {
  this->edge_container_ = g.edge_container_;
  this->expanded_edges_ = g.expanded_edges_;
  this->nodes_ = g.nodes_;
  this->id_ = g.id_;
}

ReducedGraph& ReducedGraph::operator=(const ReducedGraph& g) {
  this->expanded_edges_ = g.expanded_edges_;
  this->expanded_edges_ = g.expanded_edges_;
  this->nodes_ = g.nodes_;
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

vector<vector<Edge>> ReducedGraph::expandEdge(Edge ed){
  assert(expanded_edges_.count(ed));
  vector<vector<Edge>> all_edges;
  auto chains = expanded_edges_[ed];
  for(auto vertices : chains ){
    vector<Edge> edges;
    for(auto index=0; index<(vertices.size()-1);++index){
      Edge ed(vertices.at(index),vertices.at(index+1));
      edges.push_back(ed);
    }
    all_edges.push_back(edges);
  }
  return all_edges;
}

vector<pair<int,GraphNode>> ReducedGraph::getNodes(){
  auto vertices = edge_container_.getVertices();
  vector<pair<int,GraphNode>> nodes;
  for(auto vertex : vertices ){
    pair<int,GraphNode> pr(vertex, nodes_[vertex]);
    nodes.push_back(pr);
  }
  return nodes;
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
