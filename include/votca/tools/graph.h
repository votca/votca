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
  // First int is the index for the graph nodes, these are the same
  // indices seen in the edge container.
  std::unordered_map<int, GraphNode> nodes_;
  std::string id_;
  bool id_set_;
  void updateIds_(Graph& g);

 protected:
  void calcId_();

 public:
  Graph(){};
  /// Constructor
  Graph(std::vector<Edge> edgs, std::unordered_map<int, GraphNode> nodes)
      : EdgeContainer::EdgeContainer(edgs), nodes_(nodes), id_set_(false) {}

  std::vector<std::pair<int, GraphNode>> getIsolatedNodes(void);
  std::vector<int> getVerticesMissingNodes(void);
  std::vector<std::pair<int, GraphNode>> getNeighNodes(int vert);
  GraphNode& Node(int vert);

  GraphNode getNode(int vert);
  std::vector<std::pair<int, GraphNode>> getNodes(void);
  bool operator!=(Graph& g);
  bool operator==(Graph& g);

  Graph& operator=(const Graph& g);

  std::string getId(void);

  friend std::ostream& operator<<(std::ostream& os, const Graph g);
};

// This function is meant to be used with the stl sort algorithm e.g.:
//
// vector<pair<int,GraphNode> > vec_pr_gn = { pr_grn1, pr_grn2 , ... etc };
// sort(vec_pr_gn.begin(),vec_pr_gn.end(),cmpVertNodePairStrId);
//
bool cmpVertNodePairStrIdLessThan(std::pair<int, GraphNode> gn1_pr,
                                  std::pair<int, GraphNode> gn2_pr);
}
}
#endif  // _VOTCA_TOOLS_GRAPH_H
