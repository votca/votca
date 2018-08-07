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
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <votca/tools/edgecontainer.h>
#include <votca/tools/graphnode.h>
#include <votca/tools/name.h>

#ifndef _VOTCA_TOOLS_GRAPH_H
#define _VOTCA_TOOLS_GRAPH_H

namespace votca {
namespace tools {

/**
 * \brief A graph object that contains the graph nodes and the edges describing
 *        the bonds between nodes.
 *
 */

class GraphNode;

class Graph : public EdgeContainer {
 private:
  /// Parameter description
  /// @param int - is the index of the graph nodes / vertex ids
  /// @param GraphNode - this is the node object at each vertex and contains
  /// all the informatino that is relevant to that object
  std::unordered_map<int, GraphNode> nodes_;

  /// This is the id of the graph to graphs that contain the same content
  /// are considered equal
  std::string id_;

 protected:
  /// Calculate the id of the graph
  void calcId_();

 public:
  Graph() : id_(""){};
  ~Graph(){};

  /// Constructor
  /// @param edgs - vector of edges where each edge is composed of two
  /// ints (vertex ids) describing a link between the vertices
  /// @param nodes - unordered_map where the key is the vertex id and the
  /// target is the graph node
  Graph(std::vector<Edge> edgs, std::unordered_map<int, GraphNode> nodes)
      : EdgeContainer::EdgeContainer(edgs), nodes_(nodes) {
    calcId_();
  }

  /// Copy constructor
  Graph(const Graph& g);

  /// Equivalence and non equivalence operators work by determine if the
  /// contents of each graph node in each of the graphs are the same.
  bool operator!=(const Graph& g) const;
  bool operator==(const Graph& g) const;

  /// Copy Assignment
  Graph& operator=(const Graph& g);

  /// Move Assignment
  Graph& operator=(Graph&& g);

  /// Find all the vertices that are isolated (not connected to any other
  /// vertex) and return them in a vector with their corresponding graph node.
  std::vector<std::pair<int, GraphNode>> getIsolatedNodes(void);

  /// Functions determines which vertices do not have a graph node associated
  /// with them and return their ids in a vector.
  std::vector<int> getVerticesMissingNodes(void);

  /// Returns a vector of the vertices and their graph nodes that are directly
  /// connected to the vertex 'vert'
  std::vector<std::pair<int, GraphNode>> getNeighNodes(int vert);

  /// set the Node associated with vertex 'vert'
  void setNode(int vert, GraphNode gn);
  void setNode(std::pair<int, GraphNode> p_gn);

  /// Return a copy of the graph node at vertex 'vert'
  GraphNode getNode(int vert);

  /// Return all the vertices and their graph nodes that are within the graph
  std::vector<std::pair<int, GraphNode>> getNodes(void);

  std::string getId() { return id_; }

  friend std::ostream& operator<<(std::ostream& os, const Graph g);
};

/**
 * \brief Compare function pair<int,GraphNode> object
 *
 * This function is meant to be used with the stl sort algorithm. It will sort
 * a vector of pairs containing the vertex ids and the graphnodes. Only the
 * contetns of the graph node object are used to determine precidence e.g.
 *
 * pair<int,GraphNode> pr_grn1{ 1, gn };
 * pair<int,GraphNode> pr_grn2{ 2, gn2 };
 *
 * vector<pair<int,GraphNode> > vec_pr_gn = { pr_grn1, pr_grn2 , ... etc };
 *
 * sort(vec_pr_gn.begin(),vec_pr_gn.end(),cmpVertNodePair);
 */
bool cmpVertNodePair(std::pair<int, GraphNode> gn1_pr,
                     std::pair<int, GraphNode> gn2_pr);
}
}
#endif  // _VOTCA_TOOLS_GRAPH_H
