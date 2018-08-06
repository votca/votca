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

/**
 * \brief This file is a compilation of graph related algorithms
 *
 * These algorithms require the interplay of the graph and graph visitor
 * classes and thus cannot be made methods of either. In most cases a graph
 * visitor is specified explores the graph to determine some metric
 */
namespace votca {
namespace tools {

class Graph;
class GraphVisitor;

/// The purpose of this algorithm is to simply determine if the graph is one
/// large network or not. If it is it means every vertex is connected to 
/// every other one by 1 or more series of edges.
bool singleNetwork(Graph g, GraphVisitor& gv);

std::vector<shared_ptr<Graph>> decoupleIsolatedGraphs(Graph g);

/// This function will simply explore a graph, any information gained from the
/// exploration will depend on the graph visitor used. Note that the
/// Graph visitor is the base class which will not work on its own. The purpose
/// of doing this is to make use of polymorphism.  
void exploreGraph(Graph& g, GraphVisitor& gv);

/// This algorithm is designed to explore the topology of the graph and
/// return an identifier in the form of the string that is unique to the
/// topology. 
template <typename GV>
std::string findStructureId(Graph& g) {

  // Determine the highest degree in the graph
  int maxD = g.getMaxDegree();
  // Get the vertices with this degree
  auto verts = g.getVerticesDegree(maxD);

  // Get the nodes and determine which node has the greatest stringID
  // When compared using compare function
  std::string str_id = "";
  std::vector<int> gn_ids;
  for (auto v : verts) {
    auto gn = g.getNode(v);
    int comp_int = str_id.compare(gn.getStringId());
    if (comp_int > 0) {
      str_id = gn.getStringId();
      gn_ids.clear();
      gn_ids.push_back(v);
    } else if (comp_int == 0) {
      gn_ids.push_back(v);
    }
  }

  // If the str_id is empty it means the nodes are empty and we will
  // simply have to rely on the degree to choose the vertices to explore from
  if (str_id.compare("") == 0) {
    gn_ids = verts;
  }

  // If two or more graph nodes are found to be equal then
  // they must all be explored
  std::string chosenId = "";
  Graph g_chosen = g;

  for (auto v : gn_ids) {
    GV gv;
    gv.setStartingVertex(v);
    Graph g_temp = g;
    exploreGraph(g_temp, gv);
    std::string temp_struct_id = g_temp.getId();
    if (chosenId.compare(temp_struct_id) < 0) {
      chosenId = temp_struct_id;
      g_chosen = g_temp;
    }
  }
  
  g = g_chosen;
  return chosenId;
}
}
}

#endif  // __VOTCA_TOOLS_GRAPH_ALGORITHMS_H
