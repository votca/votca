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

#ifndef __VOTCA_TOOLS_GRAPH_VISITOR_H
#define __VOTCA_TOOLS_GRAPH_VISITOR_H

#include <iostream>
#include <set>
#include <vector>
#include <votca/tools/edge.h>
#include <votca/tools/graphnode.h>

/**
 * \brief A graph visitor class for creating graph visitor objects
 *
 * This class serves as the base class for creating graph visitor objects
 * it contains the framework by which all derived graph visitors should follow.
 * The function of the graph visitor is to explore a graph, this can be done in
 * different ways such as a breadth first or depth first algorithm or others.
 *
 * The general flow of the class is as follows:
 *
 * 1. Set the node the graph visitor will start on and create an edge queue
 * 2. Get an edge from the edge queue
 * 3. Explore the edge and operate on the unexplored node(optional) add more
 * edges to the edge queue.
 * 4. Repeat until no more edges exist in the queue.
 */
namespace votca {
namespace tools {

class Graph;

class GraphVisitor {
 protected:
  /// set containing all the vertix ids that have been explored
  std::set<int> explored_;

  /// The vertex the visitor started on
  int startingVertex_ = 0;

  /// What is done to an individual graph node as it is explored
  virtual void addEdges_(const Graph& graph, int vertex);
  virtual Edge getEdge_(const Graph& graph);
  /// Edge(0,0) is a dummy value
 public:
  virtual void exploreNode(std::pair<int, GraphNode>& vertex_and_node,
                           Graph& graph, Edge edge = DUMMY_EDGE);

  GraphVisitor(){};

  /// Determine which vertices in the edge, if any, have not been explored
  std::vector<int> getUnexploredVertex(const Edge edge) const;

  /// Determine if the exploration is complete, this is determined by whether
  /// the edge queue is empty or not, it does not necessarily mean all
  /// vertices in a graph have been explored.
  virtual bool queEmpty() const;

  void setStartingVertex(int vertex) { startingVertex_ = vertex; }
  int getStartingVertex() const { return startingVertex_; }

  /// Initialize the graphvisitor the default starting point is 0
  void initialize(Graph& graph);

  /// What the visitor does to each node as it is visited, it will
  /// simply add the vertex that was explored to the list of explored
  /// vertices in its current form.
  virtual void exec(Graph& graph, Edge edge);

  /// The next edge to be explored, note that when this function
  /// is called it removes the edge from the visitors queue and will
  /// no longer be accessible with a second call to nextEdge
  Edge nextEdge(Graph graph);

  /// Get the set of all the vertices that have been explored
  std::set<int> getExploredVertices() const;

  /// Has the vertex been explored
  bool vertexExplored(const int vert) const;
};
}  // namespace tools
}  // namespace votca

#endif  // __VOTCA_TOOLS_GRAPH_VISITOR_H
