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

#include <string>
// List of all the graph algorithms 

namespace votca {
namespace tools {

class Graph;
class GraphBasicVisitor;

bool singleNetwork(Graph g, GraphBasicVisitor& gv_b);

void exploreGraph(Graph& g, GraphBasicVisitor& gv_b, int starting_vertex=0);

std::string findStructureId(Graph& g, GraphBasicVisitor gv_b);

}
}
#endif // __VOTCA_TOOLS_GRAPH_ALGORITHMS_H
