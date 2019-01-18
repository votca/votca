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

#include <votca/csg/beadmotif.h>

#include <votca/tools/graphalgorithm.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

/**********************
 * Internal Functions *
 **********************/
bool BeadMotif::junctionExist_() {
  if (!junctionsUpToDate_) {
    junctions_ = graph_.getJunctions();
    junctionsUpToDate_ = true;
  }
  return junctions_.size() != 0;
}

bool BeadMotif::isSingle_() {
  if (BeadCount() == 1) return true;
  return false;
}

bool BeadMotif::isLine_() {
  if (junctionExist_()) return false;

  // Ensure that the degree of two of the vertices is 1
  // all other vertices must be 2
  int num_vertices_degree_1 = 0;

  vector<int> vertices = reduced_graph_.getVertices();
  for (int& vertex : vertices) {
    int degree = reduced_graph_.getDegree(vertex);
    if (degree == 1) {
      ++num_vertices_degree_1;
    } else if (degree == 0) {
      return false;
    } else if (degree > 2) {
      return false;
    }
  }
  if (num_vertices_degree_1 != 2) return false;
  return true;
}

bool BeadMotif::isLoop_() {
  if (junctionExist_()) return false;

  // Ensure that the degree of every vertex is 2
  vector<int> vertices = graph_.getVertices();
  for (int& vertex : vertices) {
    if (graph_.getDegree(vertex) != 2) return false;
  }
  return true;
}

bool BeadMotif::isFusedRing_() {
  if (!junctionExist_()) return false;
  // Ensure that the degree of every vertex is 2 or greater
  vector<int> vertices = graph_.getVertices();
  for (int& vertex : vertices) {
    if (graph_.getDegree(vertex) < 2) return false;
  }
  // If exploring from one branch of a junction lets you explore every
  // edge than it is a fused ring.
  junctionExist_();

  for (int junction : junctions_) {
    vector<Edge> edges = reduced_graph_.getNeighEdges(junction);
    set<Edge> all_edges_explored =
        exploreBranch(reduced_graph_, junction, edges.at(0));
    for (const Edge& edge_next_to_junction : edges) {
      if (!all_edges_explored.count(edge_next_to_junction)) {
        return false;
      }
    }
  }
  return true;
}
/***************************
 * Public Facing Functions *
 ***************************/

BeadMotif::MotifType BeadMotif::getType() {
  if (!type_up_to_date_) {
    CalculateType_();
  }
  return type_;
}

void BeadMotif::CalculateType_() {
  if (BeadCount() == 0) {
    type_ = MotifType::empty;
  } else if (isSingle_()) {
    type_ = MotifType::single_bead;
  } else if (!BeadStructure::isSingleStructure()) {
    type_ = MotifType::multiple_structures;
  } else if (isLine_()) {
    type_ = MotifType::line;
  } else if (isLoop_()) {
    type_ = MotifType::loop;
  } else if (isFusedRing_()) {
    type_ = MotifType::fused_ring;
  } else {
    type_ = MotifType::single_structure;
  }
  type_up_to_date_ = true;
}

BaseBead* BeadMotif::getBead(int id) { return BeadStructure::getBead(id); }

void BeadMotif::AddBead(BaseBead* bead) {
  type_ = MotifType::undefined;
  BeadStructure::AddBead(bead);
  junctionsUpToDate_ = false;
  type_up_to_date_ = false;
}

bool BeadMotif::isStructureEquivalent(BeadMotif& beadmotif) {
  return BeadStructure::isStructureEquivalent(beadmotif);
}

std::vector<BaseBead*> BeadMotif::getNeighBeads(int index) {
  return getNeighBeads(index);
}

void BeadMotif::InitializeGraph_() {
  BeadStructure::InitializeGraph_();
  reduced_graph_ = reduceGraph(graph_);
}

bool BeadMotif::isMotifSimple() {
  MotifType motif_type = getType();
  if (motif_type == single_structure || motif_type == multiple_structures ||
      motif_type == undefined) {
    return false;
  }
  return true;
}

int BeadMotif::BeadCount() { return BeadStructure::BeadCount(); }

void BeadMotif::ConnectBeads(int bead1_id, int bead2_id) {
  BeadStructure::ConnectBeads(bead1_id, bead2_id);
  junctionsUpToDate_ = false;
  type_up_to_date_ = false;
}
}  // namespace csg
}  // namespace votca
