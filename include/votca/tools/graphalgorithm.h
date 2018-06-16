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

#ifndef __VOTCA_TOOLS_GRAPH_ALGORITHMS_H
#define __VOTCA_TOOLS_GRAPH_ALGORITHMS_H

#include <iostream>
#include <string>
#include <votca/tools/graphnode.h>
// List of all the graph algorithms 

namespace votca {
namespace tools {

class Graph;
class GraphVisitor;

bool singleNetwork(Graph g, GraphVisitor& gv);

void exploreGraph(Graph& g, GraphVisitor& gv, int starting_vertex=0);

template<typename GV> 
std::string findStructureId(Graph& g){

  // Determine the highest degree in the graph
  int maxD = g.getMaxDegree();
  // Get the vertices with this degree
  auto verts = g.getVerticesDegree(maxD);

  // Get the nodes and determine which node has the greatest stringID
  // When compared using compare function
  std::string str_id = "";
  std::vector<int> gn_ids;
  for( auto v : verts ){
    auto gn = g.getNode(v);
    int comp_int = str_id.compare(gn.getStringId());
    if(comp_int>0){
      str_id = gn.getStringId();
      gn_ids.clear();
      gn_ids.push_back(v);
    }else if(comp_int==0){
      gn_ids.push_back(v);
    } 
  }

  // If the str_id is empty it means the nodes are empty and we will
  // simply have to rely on the degree to choose the vertices to explore from
  if(str_id.compare("")==0){
    gn_ids = verts;
  }

  // If two or more graph nodes are found to be equal then
  // they must all be explored   
  std::string chosenId = "";
  Graph g_chosen = g;

  for(auto v : gn_ids ){

    GV gv;
    Graph g_temp = g;
    exploreGraph(g_temp,gv,v);
    std::string temp_struct_id = g_temp.getId();
    if(chosenId.compare(temp_struct_id)>0){
      chosenId = temp_struct_id;
      g_chosen = g_temp;      
    }    
  } 

  g = g_chosen;
  return chosenId;
}

}
}
#endif // __VOTCA_TOOLS_GRAPH_ALGORITHMS_H
