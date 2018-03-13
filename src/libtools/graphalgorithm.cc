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
#include <votca/tools/graphvisitor.h>
#include <votca/tools/graphbasicvisitor.h>
#include <votca/tools/graph.h>
#include <votca/tools/graphalgorithm.h>
// List of all the graph algorithms 

using namespace std;

namespace votca {
namespace tools {

bool singleNetwork(Graph g, GraphBasicVisitor& gv_b){

  exploreGraph(g,gv_b);
  auto all_vertices = g.getVertices();
  auto exp_vertices = gv_b.getExploredVertices(); 
  auto iso_nodes = g.getIsolatedNodes();
  return exp_vertices.size()==all_vertices.size() && 
         iso_nodes.size() == 0;
}

void exploreGraph(Graph& g, GraphBasicVisitor& gv_b, int starting_vertex){
  // Create a list of all vertices and determine if they have all been
  // explored when the visitor que is empty 
  gv_b.startingVertex(g,starting_vertex);

  while(gv_b.queEmpty()==false){
    auto ed = gv_b.nextEdge(g);
    gv_b.exec(g,ed);
  }
}
/*
std::string findStructureId(Graph& g, GraphVisitor& gv_b){

  // Determine the node with the largest number of degrees
  //auto nodes = g.getNodes();
    
  // from the graph nodes sort them so that we start with the most
  // unique one first 
  // gv.startingNode(starting_node_);

    
  gv.exec();
  Edge ed = gv.nextEdge(g);

}*/

}
}
