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

#include <iostream>
#include <exception>
#include <vector>
#include <votca/tools/graphvisitor.h>
#include <votca/tools/edge.h>
#include <votca/tools/graph.h>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;
  
bool GraphVisitor::queEmpty(){
  return true;
}

void GraphVisitor::addEdges_(Graph& g, int vertex){
  throw runtime_error("addEdges_ method must be defined by your visitor");
}

void GraphVisitor::exploreNode_(pair<int,GraphNode> &p_gn,Graph& g,Edge ed){
  explored_.insert(p_gn.first);
}

vector<int> GraphVisitor::getUnexploredVertex_(Edge ed){
  vector<int> unexp_vert;
  if(explored_.count(ed.getEndPoint1())==0){
    unexp_vert.push_back(ed.getEndPoint1());
  }
  if(explored_.count(ed.getEndPoint2())==0){
    unexp_vert.push_back(ed.getEndPoint2());
  }
  return unexp_vert;
}

void GraphVisitor::initialize(Graph& g){
  auto neigh_eds = g.getNeighEdges(startingVertex_);
  GraphNode gn = g.getNode(startingVertex_);
  pair<int, GraphNode> p_gn(startingVertex_,gn);
  exploreNode_(p_gn,g);
  addEdges_(g, startingVertex_);
}

void GraphVisitor::exec(Graph& g, Edge ed){
  auto unexp_vert = getUnexploredVertex_(ed);    
  // If no vertices are return than just ignore it means the same
  // vertex was explored from a different direction
  if(!unexp_vert.size()) return;
  // If two values are returned this is a problem 
  if(unexp_vert.size()>1){
    throw runtime_error("More than one unexplored vertex in an edge,"
      " did you set the starting node");
  }
  pair<int,GraphNode> pr(unexp_vert.at(0),g.getNode(unexp_vert.at(0)));
 
  exploreNode_(pr,g,ed);

}

Edge GraphVisitor::getEdge_(Graph g){
  throw runtime_error("The getEdge_ function must be set");
}

Edge GraphVisitor::nextEdge(Graph g){

  // Get the edge and at the same time remove it from whatever queue it is in
  Edge ed = getEdge_(g);
  auto vert_v = getUnexploredVertex_(ed);
  // Do not add neighboring edges if they belong to a vertex that has already 
  // been explored because they will have already been added
  if(vert_v.size()){
    addEdges_(g, vert_v.at(0));  
  }
  return ed;
}

set<int> GraphVisitor::getExploredVertices(){ return explored_; }

}
}
