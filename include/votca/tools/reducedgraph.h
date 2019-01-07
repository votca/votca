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

class ReducedGraph : public Graph {
 private:

   std::unordered_map<Edge,std::vector<std::vector<int>>> expanded_edges_;  


 public:
  ReducedGraph(){};
  ~ReducedGraph(){};

  ReducedGraph(std::vector<ReducedEdge> reduced_edges);

  ReducedGraph(const ReducedGraph & reduced_graph);

  ReducedGraph& operator=(const ReducedGraph& reduced_graph);

  ReducedGraph& operator=(ReducedGraph&& reduced_graph);

  std::vector<Edge> getEdges();

  friend std::ostream& operator<<(std::ostream& os, const ReducedGraph g);
};

}
}
#endif  // _VOTCA_TOOLS_REDUCEDGRAPH_H
