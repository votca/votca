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

#include <votca/tools/edge.h>
#include <votca/tools/graph.h>
#include <votca/tools/graph_df_visitor.h>
#include <votca/tools/graphnode.h>

using namespace std;

namespace votca {
namespace tools {

bool Graph_DF_Visitor::queEmpty() { return edge_stack_.empty(); }

Edge Graph_DF_Visitor::getEdge_(Graph g) {
  cout << "Edge stack size " << edge_stack_.size() << endl;
  Edge ed = edge_stack_.top();
  edge_stack_.pop();
  return ed;
}

// Add edges to be explored
void Graph_DF_Visitor::addEdges_(Graph& g, int vertex) {
  auto eds = g.getNeighEdges(vertex);
  cout << "Adding Edges " << endl;
  if(edge_stack_.empty()){
  // If first edges to be added
    for(auto ed : eds){
      int neigh_vert = ed.getOtherEndPoint(vertex);
      if(explored_.count(neigh_vert)==0){
        cout << ed << endl;
        edge_stack_.push(ed);
      }
    }
  }else{
    for(auto ed : eds){
      int neigh_vert = ed.getOtherEndPoint(vertex);
      if(explored_.count(neigh_vert)==0){
        cout << ed << endl;
        edge_stack_.push(ed);
      }
    }
  }
}
}
}
