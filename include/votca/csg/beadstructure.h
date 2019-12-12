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

#pragma once
#ifndef VOTCA_CSG_BEADSTRUCTURE_H
#define VOTCA_CSG_BEADSTRUCTURE_H

#include <unordered_map>
#include <votca/tools/graph.h>

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
class BeadStructure {
 public:
  ~BeadStructure() = default;

  /**
   * @brief Given indices and edges that exist are a subset of beadstructure,
   * return the sub-beadstructure
   *
   * If the edges and and indices passed in do not exist within the current
   * BeadStructure a runtime exception is thrown, as the the sub structure is
   * not actually a sub structure of the current BeadStructure.
   *
   * @param Indices of the substructure, must exist within the current
   * BeadStructure
   * @param Edges of the substructure, must exist within the current
   * BeadStructure
   *
   * @return BeadStructure which is a substructure of BeadStructure
   */
  BeadStructure getSubStructure(const std::vector<Index> &,
                                const std::vector<tools::Edge> &) const;

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
  size_t BeadCount() const noexcept { return beads_.size(); }

  /**
   * \brief add a bead to the bead structure
   *
   * The same bead cannot be added twice.
   **/
  template <class T>
  void AddBead(const T &bead);

  /**
   * \brief Create a connection between two beads in the structure
   *
   * A bead cannot be connected to itself. It also may not be connected to a
   * bead that has not yet been added to the structure.
   **/
  virtual void ConnectBeads(const Index &bead1_id, const Index &bead2_id);

  /**
   * \brief Return a vector of all the ids of the beads neighboring the index
   **/
  std::vector<Index> getNeighBeadIds(const Index &index);

  tools::Graph getGraph();

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
  bool isStructureEquivalent(BeadStructure &beadstructure);

  /// Determine if a bead exists in the structure
  bool BeadExist(Index bead_id) const { return beads_.count(bead_id); }

  /**
   * @brief Return the ids of the beads that are in the structure
   *
   * @return vector of the ids
   */
  std::vector<Index> getBeadIds() const;

 protected:
  /// Small structure to help store bead info relevant to the structure
  struct BeadInfo {
    double mass = 0.0;
    std::string name = "";
  };
  /// If extra initialization is needed by a child class when a bead is added
  /// this method can be overridden.
  virtual void UpdateOnBeadAddition_(){};
  void InitializeGraph_();
  void CalculateStructure_();
  tools::GraphNode BeadInfoToGraphNode_(const BeadInfo &);

  bool structureIdUpToDate = false;
  bool graphUpToDate = false;
  bool single_structureUpToDate_ = false;
  bool single_structure_ = false;
  std::string structure_id_ = "";
  tools::Graph graph_;
  std::set<tools::Edge> connections_;
  std::unordered_map<Index, BeadInfo> beads_;
  std::unordered_map<Index, tools::GraphNode> graphnodes_;
};

template <class T>
void BeadStructure::AddBead(const T &bead) {
  if (beads_.count(bead.getId())) {
    std::string err = "Cannot add bead with Id ";
    err += std::to_string(bead.getId());
    err += " because each bead must have a unique Id and a bead with that Id ";
    err += "already exists within the beadstructure";
    throw std::invalid_argument(err);
  }
  UpdateOnBeadAddition_();
  size_t numberOfBeads = beads_.size();
  BeadInfo bead_info;
  bead_info.mass = bead.getMass();
  bead_info.name = bead.getName();

  beads_[bead.getId()] = bead_info;

  if (numberOfBeads != beads_.size()) {
    single_structureUpToDate_ = false;
    graphUpToDate = false;
    structureIdUpToDate = false;
  }
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADSTRUCTURE_H
