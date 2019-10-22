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

#ifndef _VOTCA_CSG_BEADSTRUCTURE_H
#define _VOTCA_CSG_BEADSTRUCTURE_H

#include <cassert>
#include <unordered_map>
#include <votca/tools/graph.h>
#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graphalgorithm.h>
#include <votca/tools/graphdistvisitor.h>

namespace TOOLS = votca::tools;

namespace votca {
namespace csg {

/**
 * \brief Designed to determine if the structure beads passed in
 *
 * Essentially it will have the functionality to determine if the stored beads
 * make up a single molecule. It can also break the stored beads up into
 * molecules. It can compare two bead structures and determine if they are the
 * same structure. The Ids of a bead will not determine if a structure is
 * equivalent or not. Each bead must have a unique Id.
 *
 * E.g. If a beadstructure is composed of the following:
 *
 * BeadStructure 1
 * Name:"A" Id:1 ---- Name:"B" Id:2
 *
 * BeadStucture 2
 * Name:"A" Id:3 ---- Name:"B" Id:4
 *
 * Here BeadStructure 1 and 2 will compare equal
 *
 **/
template <class T>
class BeadStructure {
 public:
  ~BeadStructure() = default;

  /**
   * \brief Determine if the bead structure consists of a single connected
   * structure
   *
   * This function will determine if all beads in the structure are connected
   * somehow to everyother bead. The connection does not have to be direct
   *
   * @return - returns a boolean true if it is a single Structure
   **/
  bool isSingleStructure();

  /**
   * \brief returns the number of beads in the bead structure
   **/
  size_t BeadCount() { return beads_.size(); }

  /**
   * \brief add a bead to the bead structure
   *
   * The same bead cannot be added twice.
   **/
  void AddBead(T *bead);

  /**
   * \brief Get the bead with the specified id
   **/
  T *getBead(int id);

  /**
   * \brief Create a connection between two beads in the structure
   *
   * A bead cannot be connected to itself. It also may not be connected to a
   * bead that has not yet been added to the structure.
   **/
  void ConnectBeads(int bead1_id, int bead2_id);

  /**
   * \brief Return a vector of all the beads neighboring the index
   **/
  std::vector<T *> getNeighBeads(int index);

  TOOLS::Graph getGraph();
  /**
   * \brief Compare the topology of two bead structures
   *
   * This function looks at how the beads are arranged within the bead structure
   * and determines if the topology is the same.
   *
   * @param[in] - beadstructure to compare with
   * @return - if the same returns true else false
   *
   **/
  bool isStructureEquivalent(BeadStructure<T> &beadstructure);

  bool BeadExist(int bead_id) const { return beads_.count(bead_id); }

 protected:
  void InitializeGraph_();
  void CalculateStructure_();
  TOOLS::GraphNode BaseBeadToGraphNode_(T *basebead);

  bool structureIdUpToDate = false;
  bool graphUpToDate = false;
  bool single_structureUpToDate_ = false;
  bool single_structure_ = false;
  std::string structure_id_ = "";
  TOOLS::Graph graph_;
  std::set<TOOLS::Edge> connections_;
  std::unordered_map<int, T *> beads_;
  std::unordered_map<long int, TOOLS::GraphNode> graphnodes_;
};

/**********************
 * Internal Functions *
 **********************/

template <class T>
void BeadStructure<T>::InitializeGraph_() {
  if (!graphUpToDate) {
    std::vector<TOOLS::Edge> connections_vector;
    for (const TOOLS::Edge &edge : connections_) {
      connections_vector.push_back(edge);
    }

    for (std::pair<const int, T *> &id_bead_ptr_pair : beads_) {
      graphnodes_[id_bead_ptr_pair.first] =
          BaseBeadToGraphNode_(id_bead_ptr_pair.second);
    }
    graph_ = TOOLS::Graph(connections_vector, graphnodes_);
    graphUpToDate = true;
  }
}

template <class T>
void BeadStructure<T>::CalculateStructure_() {

  InitializeGraph_();
  if (!structureIdUpToDate) {
    structure_id_ = TOOLS::findStructureId<TOOLS::GraphDistVisitor>(graph_);
    structureIdUpToDate = true;
  }
}

template <class T>
TOOLS::GraphNode BeadStructure<T>::BaseBeadToGraphNode_(T *basebead) {
  std::unordered_map<std::string, double> attributes1;
  std::unordered_map<std::string, std::string> attributes2;

  attributes1["Mass"] = basebead->getMass();
  attributes2["Name"] = basebead->getName();

  /// Add graphnodes
  TOOLS::GraphNode graphnode;
  graphnode.setDouble(attributes1);
  graphnode.setStr(attributes2);

  return graphnode;
}

/***************************
 * Public Facing Functions *
 ***************************/

template <class T>
void BeadStructure<T>::AddBead(T *bead) {
  if (beads_.count(bead->getId())) {
    std::string err = "Cannot add bead with Id ";
    err += std::to_string(bead->getId());
    err += " because each bead must have a unique Id and a bead with that Id ";
    err += "already exists within the beadstructure";
    throw std::invalid_argument(err);
  }
  size_t numberOfBeads = beads_.size();
  beads_[bead->getId()] = bead;

  if (numberOfBeads != beads_.size()) {
    single_structureUpToDate_ = false;
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

template <class T>
void BeadStructure<T>::ConnectBeads(int bead1_id, int bead2_id) {
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
  connections_.insert(TOOLS::Edge(bead1_id, bead2_id));
  if (numberOfConnections != connections_.size()) {
    single_structureUpToDate_ = false;
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

template <class T>
TOOLS::Graph BeadStructure<T>::getGraph() {
  InitializeGraph_();
  return graph_;
}

template <class T>
bool BeadStructure<T>::isSingleStructure() {

  InitializeGraph_();
  if (single_structureUpToDate_ == false) {
    std::vector<long int> vertices = graph_.getVertices();
    if (vertices.size() == 0) {
      single_structure_ = false;
      return single_structure_;
    }
    // Choose first vertex that is actually in the graph as the starting vertex
    TOOLS::Graph_BF_Visitor gv_breadth_first;
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

template <class T>
bool BeadStructure<T>::isStructureEquivalent(BeadStructure<T> &beadstructure) {
  if (!structureIdUpToDate) {
    CalculateStructure_();
  }
  if (!beadstructure.structureIdUpToDate) {
    beadstructure.CalculateStructure_();
  }
  return structure_id_.compare(beadstructure.structure_id_) == 0;
}

template <class T>
std::vector<T *> BeadStructure<T>::getNeighBeads(int index) {
  if (!graphUpToDate) {
    InitializeGraph_();
  }
  std::vector<int> neighbor_ids = graph_.getNeighVertices(index);
  std::vector<T *> neighbeads;
  for (int &node_id : neighbor_ids) {
    neighbeads.push_back(beads_[node_id]);
  }
  return neighbeads;
}

template <class T>
T *BeadStructure<T>::getBead(int index) {
  assert(beads_.count(index));
  return beads_[index];
}

}  // namespace csg
}  // namespace votca

#endif
