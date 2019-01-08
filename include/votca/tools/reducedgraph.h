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

#include <votca/tools/graph.h>
#include <votca/tools/reducededge.h>

#ifndef _VOTCA_TOOLS_REDUCEDGRAPH_H
#define _VOTCA_TOOLS_REDUCEDGRAPH_H

namespace votca {
namespace tools {

/**
 * \brief Contains a graph that consits of vertices with degree of 1 or greater than 3
 *
 * The point of this class is to reduce the computational complexity of a
 * regular graph. This is achieved by removing any vertices with degree 2. For 
 * example a graph:
 *
 * 1 - 2 - 3 - 4 - 5 - 9
 *     |   |   |
 *     6 - 7   8
 *
 * Would be reduced to 
 *
 * 1 - 2 - 3 - 4 - 9
 *     | _ |   |
 *             8
 *
 * Notice that the vertices 5, 6 and 7 have been removed, there also exist two
 * edges connecting 2 to 3. 
 *
 **/
class ReducedGraph : public Graph {
 private:

   std::unordered_map<Edge,std::vector<std::vector<int>>> expanded_edges_;  
   std::unordered_map<int,GraphNode> hidden_nodes_;

 public:
  ReducedGraph(){};
  ~ReducedGraph(){};

  ReducedGraph(std::vector<ReducedEdge> reduced_edges);
  ReducedGraph(std::vector<ReducedEdge> reduced_edges, std::unordered_map<int,GraphNode> nodes);

  ReducedGraph(const ReducedGraph & reduced_graph);

  ReducedGraph& operator=(const ReducedGraph& reduced_graph);

  ReducedGraph& operator=(ReducedGraph&& reduced_graph);

  std::vector<Edge> getEdges();

  std::vector<std::vector<Edge>> expandEdge(Edge ed);

  int getDegree(int vertex);

  friend std::ostream& operator<<(std::ostream& os, const ReducedGraph g);
};

}
}
#endif  // _VOTCA_TOOLS_REDUCEDGRAPH_H
