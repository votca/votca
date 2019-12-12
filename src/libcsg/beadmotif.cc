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

#include "../../include/votca/csg/beadmotif.h"
#include <votca/tools/graphalgorithm.h>
namespace votca {
namespace csg {
class BaseBead;
}  // namespace csg
}  // namespace votca

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

/**********************
 * Internal Functions *
 **********************/

void BeadMotif::InitializeGraph_() {
  BeadStructure::InitializeGraph_();
  reduced_graph_ = reduceGraph(graph_);
}

bool BeadMotif::junctionExist_() {
  if (!junctionsUpToDate_) {
    junctions_ = graph_.getJunctions();
    junctionsUpToDate_ = true;
  }
  return junctions_.size() != 0;
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

bool BeadMotif::isSingle_() const noexcept {
  if (BeadCount() == 1) {
    return true;
  }
  return false;
}

bool BeadMotif::isLine_() {
  if (junctionExist_()) {
    return false;
  }
  // Ensure that the degree of two of the vertices is 1
  // all other vertices must be 2
  Index num_vertices_degree_1 = 0;

  vector<Index> vertices = reduced_graph_.getVertices();
  for (Index& vertex : vertices) {
    Index degree = reduced_graph_.getDegree(vertex);
    if (degree == 1) {
      ++num_vertices_degree_1;
    } else if (degree == 0) {
      return false;
    } else if (degree > 2) {
      return false;
    }
  }
  if (num_vertices_degree_1 != 2) {
    return false;
  }
  return true;
}

bool BeadMotif::isLoop_() {
  if (junctionExist_()) {
    return false;
  }

  // Ensure that the degree of every vertex is 2
  vector<Index> vertices = graph_.getVertices();
  for (Index& vertex : vertices) {
    if (graph_.getDegree(vertex) != 2) {
      return false;
    }
  }
  return true;
}

/**
 * One has to explore the whole tree from each of the junctions to
 * determine if the model is a fused ring or not. For speed it might
 * make since to reduce the graph first to junctions of 3 or more.
 *
 * if There is not way back to the junction than you have something
 * like this:
 *
 * c1 - c2 - c5 - c6
 * |    |    |    |
 * c3 - c4   c7 - c8
 *
 * Say you start at c2 and head down tree c5 if you never find a way back
 * you can split it
 *
 * If you have something like this
 *
 * c1 - c2 - c3
 * |  /   \  |
 *  c4     c5
 *
 *  Then you do not have a fused ring, must be represented as a joint
 *  and two lines. Exploring a tree will only lead to one way back.
 *
 *         c6
 *        /  |
 * c1 - c2 - c3
 * |  /   \  |
 *  c4     c5
 *
 *  Still acts like a joint, For it not to be a joint exploring a single
 *  branch originating from the junction should lead to exploration of
 *  all the edges.
 **/
bool BeadMotif::isFusedRing_() {
  if (!junctionExist_()) {
    return false;
  }
  // Ensure that the degree of every vertex is 2 or greater
  vector<Index> vertices = graph_.getVertices();
  for (Index& vertex : vertices) {
    if (graph_.getDegree(vertex) < 2) {
      return false;
    }
  }
  // If exploring from one branch of a junction lets you explore every
  // edge than it is a fused ring.
  junctionExist_();

  for (Index junction : junctions_) {
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

void BeadMotif::UpdateOnBeadAddition_() {
  type_ = MotifType::undefined;
  junctionsUpToDate_ = false;
  type_up_to_date_ = false;
}

/***************************
 * Public Facing Functions *
 ***************************/

BeadMotif::MotifType BeadMotif::getType() {
  if (!type_up_to_date_) {
    InitializeGraph_();
    CalculateType_();
  }
  return type_;
}

bool BeadMotif::isMotifSimple() {
  MotifType motif_type = getType();
  if (motif_type == single_structure || motif_type == multiple_structures ||
      motif_type == undefined) {
    return false;
  }
  return true;
}

void BeadMotif::ConnectBeads(const Index& bead1_id, const Index& bead2_id) {
  BeadStructure::ConnectBeads(bead1_id, bead2_id);
  junctionsUpToDate_ = false;
  type_up_to_date_ = false;
}

}  // namespace csg
}  // namespace votca
