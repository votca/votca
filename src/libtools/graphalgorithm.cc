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

void exploreGraph(Graph& g, GraphVisitor& gv){

  if(!g.vertexExist(gv.getStartingVertex())){
    string err = "Cannot explore graph starting at vertex " +
      to_string(gv.getStartingVertex()) + " because it does not exist in the "
      "graph. You can change the starting vertex by using the "
      "setStartingVertex method of the visitor instance.";
    throw invalid_argument(err);
  }
  // Create a list of all vertices and determine if they have all been
  // explored when the visitor que is empty 
  gv.initialize(g);
  while(!gv.queEmpty()){
    auto ed = gv.nextEdge(g);
    gv.exec(g,ed);
  }
}

}
}
