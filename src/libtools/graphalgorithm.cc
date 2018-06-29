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

#include <string>
#include <votca/tools/graphvisitor.h>
#include <votca/tools/graph.h>
#include <votca/tools/graphalgorithm.h>

using namespace std;

namespace votca {
namespace tools {

bool singleNetwork(Graph g, GraphVisitor& gv){
  exploreGraph(g,gv);
  return gv.getExploredVertices().size()==g.getVertices().size() && 
         g.getIsolatedNodes().size() == 0;
}

void exploreGraph(Graph& g, GraphVisitor& gv, int starting_vertex){

  // Create a list of all vertices and determine if they have all been
  // explored when the visitor que is empty 
  gv.startingVertex(g,starting_vertex);
  while(!gv.queEmpty()){
    auto ed = gv.nextEdge(g);
    gv.exec(g,ed);
  }
}

}
}
