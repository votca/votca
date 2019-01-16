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
#include <algorithm>
#include <votca/tools/edge.h>
#include <votca/tools/graph.h>
#include <votca/tools/graph_df_visitor.h>
#include <votca/tools/graphnode.h>

using namespace std;

namespace votca {
namespace tools {

bool Graph_DF_Visitor::queEmpty() { return edge_list_.empty(); }

Edge Graph_DF_Visitor::getEdge_(Graph g) {
  Edge ed = edge_list_.back();
  edge_list_.pop_back();
  return ed;
}

// Add edges to be explored
void Graph_DF_Visitor::addEdges_(Graph& g, int vertex) {
  auto eds = g.getNeighEdges(vertex);
  if(edge_list_.empty()){
  // If first edges to be added
    for(auto ed : eds){
      int neigh_vert = ed.getOtherEndPoint(vertex);
      if(explored_.count(neigh_vert)==0){
        edge_list_.push_back(ed);
      }
    }
  }else{
    for(const auto& ed : eds){
      int neigh_vert = ed.getOtherEndPoint(vertex);
      if(explored_.count(neigh_vert)==0){
        edge_list_.push_back(ed);
      }else{
        // Check if edge has already been added earlier in the queue
        // if so we wil move it to the end
        list<Edge>::iterator edge_found_iterator = find(edge_list_.begin(),edge_list_.end(),ed);
        if(edge_found_iterator!=edge_list_.end()){
          // Move the edge to the back if it was stored earlier on.
          edge_list_.splice(edge_list_.end(),edge_list_,edge_found_iterator); 
        }
      }
    }
  }
}
}
}
