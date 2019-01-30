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

#include <iostream>
#include <unordered_map>
#include <votca/csg/beadstructure.h>
#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graphalgorithm.h>
#include <votca/tools/graphdistvisitor.h>

using namespace votca::tools;
using namespace std;

namespace votca {

namespace csg {
/**********************
 * Internal Functions *
 **********************/

GraphNode BaseBeadToGraphNode(BaseBead *basebead) {
  unordered_map<string, double> attributes1;
  unordered_map<string, string> attributes2;

  attributes1["Mass"] = basebead->getMass();
  attributes2["Name"] = basebead->getName();

  /// Add graphnodes
  GraphNode graphnode;
  graphnode.setDouble(attributes1);
  graphnode.setStr(attributes2);

  return graphnode;
}

/***************************
 * Public Facing Functions *
 ***************************/

void BeadStructure::AddBead(BaseBead *bead) {
  if (beads_.count(bead->getId())) {
    string err = "Cannot add bead with Id ";
    err += to_string(bead->getId());
    err += " because each bead must have a unique Id and a bead with that Id ";
    err += "already exists within the beadstructure";
    throw invalid_argument(err);
  }
  size_t numberOfBeads = beads_.size();
  beads_[bead->getId()] = bead;

  if (numberOfBeads != beads_.size()) {
    single_structureUpToDate_ = false;
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

void BeadStructure::ConnectBeads(int bead1_id, int bead2_id) {
  if (!(beads_.count(bead1_id)) || !(beads_.count(bead2_id))) {
    string err = "Cannot connect beads in bead structure that do not exist";
    throw invalid_argument(err);
  }
  if (bead1_id == bead2_id) {
    string err = "Beads cannot be self-connected";
    throw invalid_argument(err);
  }
  size_t numberOfConnections = connections_.size();
  connections_.insert(Edge(bead1_id, bead2_id));
  if (numberOfConnections != connections_.size()) {
    single_structureUpToDate_ = false;
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

void BeadStructure::InitializeGraph_() {
  if (!graphUpToDate) {
    vector<Edge> connections_vector;
    for (const Edge &edge : connections_) {
      connections_vector.push_back(edge);
    }

    for (pair<const int, BaseBead *> &id_bead_ptr_pair : beads_) {
      graphnodes_[id_bead_ptr_pair.first] =
          BaseBeadToGraphNode(id_bead_ptr_pair.second);
    }
    graph_ = Graph(connections_vector, graphnodes_);
    graphUpToDate = true;
  }
}

Graph BeadStructure::getGraph() {
  InitializeGraph_();
  return graph_;
}

void BeadStructure::CalculateStructure_() {

  InitializeGraph_();
  if (!structureIdUpToDate) {
    structure_id_ = findStructureId<GraphDistVisitor>(graph_);
    structureIdUpToDate = true;
  }
}

bool BeadStructure::isSingleStructure() {

  InitializeGraph_();
  if (single_structureUpToDate_ == false) {
    vector<int> vertices = graph_.getVertices();
    if (vertices.size() == 0) {
      single_structure_ = false;
      return single_structure_;
    }
    // Choose first vertex that is actually in the graph as the starting vertex
    Graph_BF_Visitor gv_breadth_first;
    gv_breadth_first.setStartingVertex(vertices.at(0));
    if (!isSingleNetwork(graph_, gv_breadth_first)) {
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

vector<BaseBead *> BeadStructure::getNeighBeads(int index) {
  if (!graphUpToDate) {
    InitializeGraph_();
  }
  vector<int> neighbor_ids = graph_.getNeighVertices(index);
  vector<BaseBead *> neighbeads;
  for (int &node_id : neighbor_ids) {
    neighbeads.push_back(beads_[node_id]);
  }
  return neighbeads;
}

BaseBead *BeadStructure::getBead(int index) {
  assert(beads_.count(index));
  return beads_[index];
}

}  // namespace csg
}  // namespace votca
