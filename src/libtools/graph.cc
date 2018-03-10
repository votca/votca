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
#include <votca/tools/graph.h>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;

void Graph::updateStructIds_(Graph& g){
  if(!this->structure_id_set_){
    this->calcStructureId_();
    this->structure_id_set_=true;
  }
  if(!g.structure_id_set_){
    g.calcStructureId_();
    g.structure_id_set_=true;
  }
}

bool Graph::operator!=(Graph& g) {
  updateStructIds_(g);
  return structure_id_.compare(g.structure_id_);
}

bool Graph::operator==( Graph& g) {
  return !(*(this)!=g);  
}

GraphNode& Graph::getNode(int vert){
  return nodes_[vert];
}

void Graph::calcStructureId_(){
  throw runtime_error("calcStructureId_ is not yet implemented");
  return;
}

}
}

