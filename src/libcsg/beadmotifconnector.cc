/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "../../include/votca/csg/beadmotifconnector.h"

using namespace std;
using namespace votca::tools;

namespace votca {
namespace csg {

using element = reduced_edge_to_edges_map::value_type;

void BeadMotifConnector::AddMotifAndBeadEdge(Edge motif_edge, Edge bead_edge) {
  motif_and_bead_edges_.insert(element(motif_edge, bead_edge));
}

vector<Edge> BeadMotifConnector::getBeadEdges(Edge motif_edge) {
  auto left_iterator_range = motif_and_bead_edges_.left.equal_range(motif_edge);
  vector<Edge> bead_edges;
  for (auto iterator = left_iterator_range.first;
       iterator != left_iterator_range.second; ++iterator) {
    bead_edges.push_back(iterator->second);
  }
  return bead_edges;
}

vector<Edge> BeadMotifConnector::getBeadEdges() {
  vector<Edge> bead_edges;
  for (auto left_iterator = motif_and_bead_edges_.left.begin();
       left_iterator != motif_and_bead_edges_.left.end(); ++left_iterator) {
    bead_edges.push_back(left_iterator->second);
  }
  return bead_edges;
}

Edge BeadMotifConnector::getMotifEdge(Edge bead_edge) {
  return motif_and_bead_edges_.right.at(bead_edge);
}

unordered_set<Edge> BeadMotifConnector::getMotifEdges() {
  unordered_set<Edge> motif_edges;
  for (auto left_iterator = motif_and_bead_edges_.left.begin();
       left_iterator != motif_and_bead_edges_.left.end(); ++left_iterator) {
    motif_edges.insert(left_iterator->first);
  }
  return motif_edges;
}
}  // namespace csg
}  // namespace votca
