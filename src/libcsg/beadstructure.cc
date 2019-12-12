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
#include <cassert>

#include "../../include/votca/csg/beadstructure.h"

#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graphalgorithm.h>
#include <votca/tools/graphdistvisitor.h>

using namespace std;
using namespace votca;
using namespace votca::csg;
using namespace votca::tools;

/**********************
 * Internal Functions *
 **********************/

void BeadStructure::InitializeGraph_() {
  if (!graphUpToDate) {
    std::vector<tools::Edge> connections_vector;
    for (const tools::Edge &edge : connections_) {
      connections_vector.push_back(edge);
    }

    for (std::pair<const Index, BeadInfo> &id_bead_ptr_pair : beads_) {
      graphnodes_[id_bead_ptr_pair.first] =
          BeadInfoToGraphNode_(id_bead_ptr_pair.second);
    }
    graph_ = tools::Graph(connections_vector, graphnodes_);
    graphUpToDate = true;
  }
}

void BeadStructure::CalculateStructure_() {

  InitializeGraph_();
  if (!structureIdUpToDate) {
    structure_id_ = tools::findStructureId<tools::GraphDistVisitor>(graph_);
    structureIdUpToDate = true;
  }
}

tools::GraphNode BeadStructure::BeadInfoToGraphNode_(
    const BeadInfo &bead_info) {
  std::unordered_map<std::string, double> attributes1;
  std::unordered_map<std::string, std::string> attributes2;

  attributes1["Mass"] = bead_info.mass;
  attributes2["Name"] = bead_info.name;

  /// Add graphnodes
  tools::GraphNode graphnode;
  graphnode.setDouble(attributes1);
  graphnode.setStr(attributes2);

  return graphnode;
}

/***************************
 * Public Facing Functions *
 ***************************/

BeadStructure BeadStructure::getSubStructure(
    const std::vector<Index> &bead_ids,
    const std::vector<tools::Edge> &connections) const {
  BeadStructure new_beadstructure;
  for (const Index &bead_id : bead_ids) {
    new_beadstructure.beads_[bead_id] = beads_.at(bead_id);
  }
  for (const tools::Edge &edge : connections) {
    new_beadstructure.ConnectBeads(edge.getEndPoint1(), edge.getEndPoint2());
  }
  return new_beadstructure;
}

void BeadStructure::ConnectBeads(const Index &bead1_id, const Index &bead2_id) {
  if (!(beads_.count(bead1_id)) || !(beads_.count(bead2_id))) {
    std::string err =
        "Cannot connect beads in bead structure that do not exist";
    throw std::invalid_argument(err);
  }
  if (bead1_id == bead2_id) {
    std::string err = "Beads cannot be self-connected";
    throw std::invalid_argument(err);
  }
  size_t numberOfConnections = connections_.size();
  connections_.insert(tools::Edge(bead1_id, bead2_id));
  if (numberOfConnections != connections_.size()) {
    single_structureUpToDate_ = false;
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

tools::Graph BeadStructure::getGraph() {
  InitializeGraph_();
  return graph_;
}

bool BeadStructure::isSingleStructure() {

  InitializeGraph_();
  if (single_structureUpToDate_ == false) {
    std::vector<Index> vertices = graph_.getVertices();
    if (vertices.size() == 0) {
      single_structure_ = false;
      return single_structure_;
    }
    // Choose first vertex that is actually in the graph as the starting vertex
    tools::Graph_BF_Visitor gv_breadth_first;
    gv_breadth_first.setStartingVertex(vertices.at(0));
    if (!singleNetwork(graph_, gv_breadth_first)) {
      single_structure_ = false;
      return single_structure_;
    }
    if (beads_.size() == 0) {
      single_structure_ = false;
      return single_structure_;
    }
    if (vertices.size() != beads_.size()) {
      single_structure_ = false;
      return single_structure_;
    }
    single_structure_ = true;
    single_structureUpToDate_ = true;
  }
  return single_structure_;
}

bool BeadStructure::isStructureEquivalent(BeadStructure &beadstructure) {
  if (!structureIdUpToDate) {
    CalculateStructure_();
  }
  if (!beadstructure.structureIdUpToDate) {
    beadstructure.CalculateStructure_();
  }
  return structure_id_.compare(beadstructure.structure_id_) == 0;
}

std::vector<Index> BeadStructure::getNeighBeadIds(const Index &index) {
  if (!graphUpToDate) {
    InitializeGraph_();
  }
  std::vector<Index> neighbor_ids = graph_.getNeighVertices(index);
  return neighbor_ids;
}

std::vector<Index> BeadStructure::getBeadIds() const {
  /// Do not use graph_.getVertices() because this would require that the graph
  /// is updated
  vector<Index> bead_ids;
  for (auto &id_and_bead_info : beads_) {
    bead_ids.push_back(id_and_bead_info.first);
  }
  return bead_ids;
}
