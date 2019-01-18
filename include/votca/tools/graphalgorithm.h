/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#include <memory>
#include <string>
#include <votca/tools/graphnode.h>
#include <votca/tools/reducedgraph.h>

/**
 * \brief This file is a compilation of graph related algorithms.
 *
 * These algorithms require the interplay of the graph and graph visitor
 * classes and thus cannot be made methods of either. In most cases a graph
 * visitor is specified explores the graph to determine some metric
 */
namespace votca {
namespace tools {

class Graph;
class GraphVisitor;

/**
 * \brief Determine if every vertex is connected to every other one through some
 *        combination of edges.
 *
 * The purpose of this algorithm is to simply determine if the graph is one
 * large network or not. If it is it means every vertex is connected to every
 * other one by 1 or more series of edges.
 *
 * @param[in] - Graph instance
 * @param[in,out] - Graph visitor reference instance used to explore the graph
 * @return - Boolean value (true - if single network)
 */
bool singleNetwork(Graph graph, GraphVisitor& graph_visitor);

/**
 * \brief Explore one of the branches if exploration is initiated at vertex
 * `starting_vertex`.
 *
 * The `edge` is used to determine which branch is to be explored. The edge must
 * be an edge connected to the starting_vertex.
 *
 * Take the following
 *
 * 1 - 2 - 3 - 4 - 5
 *         |       |
 *         6 - 7 - 8
 *
 * If our graph is reprsented by the above depiction and vertex 3 is chosen as
 * our starting vertex, we are left with 3 edges to choose from Edges:
 *
 *   2 - 3
 *   3 - 4
 *   3 - 6
 *
 * If we pick edge 3 - 4 or 3 - 6 the returned set of edges will consist of the
 * loop
 *
 *         3 - 4 - 5
 *         |       |
 *         6 - 7 - 8
 *
 * If instead Edge 2 - 3 is picked, the following set of edges would be returned
 *
 * 1 - 2 - 3
 *
 * @param[in] - Graph instance
 * @param[in] - int starting vertex, where the exploration begins
 * @param[in] - the edge indicating which branch is to be explored
 * @return - set of edges in the branch that were explored
 **/
std::set<Edge> exploreBranch(Graph g, int starting_vertex, Edge edge);

/**
 * \brief Will take a graph and reduce it, by removing all vertices with degree
 * of 2.
 *
 * The purpose of this algorithm is to introduce new functionality that reduces
 * the complexity of a graph by removing any vertex that has a degree of 2. By
 * exploring a reduced graph instead of a full graph insight can still be gained
 * into the topology of the graph but with better performance. The edges of the
 * reduced graph can be expanded if needed to determine how they correspond to
 * the full graph. Take:
 *
 * 1 - 2 - 3 - 4 - 5 - 6
 *     |       |
 *     7 - 8 - 9
 *
 * This would be reduced to
 *
 * 1 - 2 - 4 - 6
 *     | _ |
 *
 * A total of 4 vertices with 4 edges as opposed to 9 vertices and 9 edges.
 *
 * @param[in] - graph instance
 * @return - a reduced graph
 **/
ReducedGraph reduceGraph(Graph graph);

/**
 * \brief Break graph into smaller graph instances if the network is made up of
 *        isolated sub networks.
 *
 * This algorithm will determine if there are groups of vertices that are
 * connected, but where there are no connections shared between the groups.
 * These groups will be brocken up into their own Graph instances.
 *
 * @param[in] - Graph instance
 * @return - vector containing shared pointers to all the sub graphs if there
 *           are no subgraphs than the input graph is returned.
 */
std::vector<Graph> decoupleIsolatedSubGraphs(Graph& graph);

/**
 * \brief Explore a graph with a graph visitor.
 *
 * This function will simply explore a graph, any information gained from the
 * exploration will depend on the graph visitor used. Note that the Graph
 * visitor is the base class which will not work on its own. The purpose of
 * doing this is to make use of polymorphism.
 *
 * @param[in,out] - Graph reference instance
 * @param[in,out] - graph visitor
 */
void exploreGraph(Graph& graph, GraphVisitor& graph_visitor);

/**
 * \brief Find a unique identifier that describes graph structure.
 *
 * This algorithm is designed to explore the topology of the graph and return an
 * identifier in the form of the string that is unique to the topology. It does
 * this by taking into account the contents of the graphnodes. How it does this
 * is specific to the graph visitor specified.
 *
 * @param[in,out] - Graph reference instance
 * @return - string identifier
 */
template <typename GV>
std::string findStructureId(Graph& graph) {

  // Determine the highest degree in the graph
  int maxD = graph.getMaxDegree();
  // Get the vertices with this degree
  std::vector<int> vertices = graph.getVerticesDegree(maxD);

  // Get the nodes and determine which node has the greatest stringID
  // When compared using compare function
  std::string str_id = "";
  std::vector<int> graph_node_ids;
  for (const int& vertex : vertices) {
    auto graph_node = graph.getNode(vertex);
    int comp_int = str_id.compare(graph_node.getStringId());
    if (comp_int > 0) {
      str_id = graph_node.getStringId();
      graph_node_ids.clear();
      graph_node_ids.push_back(vertex);
    } else if (comp_int == 0) {
      graph_node_ids.push_back(vertex);
    }
  }

  // If the str_id is empty it means the nodes are empty and we will
  // simply have to rely on the degree to choose the vertices to explore from
  if (str_id.compare("") == 0) {
    graph_node_ids = vertices;
  }
  // If two or more graph nodes are found to be equal then
  // they must all be explored
  std::string chosenId = "";
  Graph graph_chosen = graph;

  for (const int& vertex : graph_node_ids) {
    GV graph_visitor;
    graph_visitor.setStartingVertex(vertex);
    Graph graph_temp = graph;
    exploreGraph(graph_temp, graph_visitor);
    std::string temp_struct_id = graph_temp.getId();
    if (chosenId.compare(temp_struct_id) < 0) {
      chosenId = temp_struct_id;
      graph_chosen = graph_temp;
    }
  }

  graph = graph_chosen;
  return chosenId;
}
}  // namespace tools
}  // namespace votca

#endif  // __VOTCA_TOOLS_GRAPH_ALGORITHMS_H
