/*
 *            Copyright 2018-2018 The VOTCA Development Team
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

#include <unordered_map>
#include <string>
#include <iostream>
#include <utility>
#include <votca/tools/name.h>
#include <vector>
#include <votca/tools/graphnode.h>
#include <votca/tools/edgecontainer.h>

#ifndef _VOTCA_TOOLS_GRAPH_H
#define _VOTCA_TOOLS_GRAPH_H

namespace votca {
namespace tools {

/**
 * \brief A graph object that contains the graph nodes and the edges describing
 *        the bonds between nodes.
 *
 */
class Graph {
  private:
    EdgeContainer edges_;
    std::unordered_map<int,GraphNode> nodes_;
  public:
  /// Constructor
  /// Basic setters
};
}
}
#endif  // _VOTCA_TOOLS_GRAPH_H
