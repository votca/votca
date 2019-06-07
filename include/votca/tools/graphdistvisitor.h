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

#pragma once
#ifndef VOTCA_TOOLS_GRAPH_DIST_VISITOR_H
#define VOTCA_TOOLS_GRAPH_DIST_VISITOR_H

#include "graph_bf_visitor.h"
#include <deque>
#include <iostream>
#include <queue>

/**
 * \brief A graph visitor determines the graph topology
 *
 * This visitor will calculate the distance of each node from the starting node
 * it is built on top of the graph breadth first visitor. As the visitor moves
 * through the graph it adds a 'Dist' attribute to each graph node with an
 * integer value corresponding to how far it is removed from the starting node.
 *
 * E.g.
 *
 * 0 - 1 - 2 - 3
 *
 * If vertex 1 is the starting vertex than the graph node associated with
 * vertex 1 will have a distance of 0. Vertices 0 and 2 a distance of 1 and
 * vertex 3 a distnace of 2.
 */
namespace votca {
namespace tools {

class Graph;
class Edge;
class GraphNode;

class GraphDistVisitor : public Graph_BF_Visitor {

 public:
  GraphDistVisitor(){};

  /// Note the only manipulation to the BF visitor is the need to add a
  /// distance attribute to each of the graph nodes.
  void exploreNode(std::pair<int, GraphNode>& p_gn, Graph& g,
                   Edge ed = DUMMY_EDGE);
};
}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_GRAPH_DIST_VISITOR_H
