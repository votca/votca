/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writingraph, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_TOOLS_BRANCHSEQUENCELEAF_H
#define VOTCA_TOOLS_BRANCHSEQUENCELEAF_H
#pragma once
#include "branch.h"
#include "votca/tools/contentlabel.h"
#include <iostream>
#include <memory>
#include <string>

namespace votca {
namespace tools {
/**
 * @brief Private class to work with canonize algorithm
 *
 * This class is designed to sort branches of a graph into a unique order.
 * The idea is that no matter the order if a unique order of the vertices
 * exists it can be found using this class. The class only relies on the
 * contents of the leafs and thier relationship to one another, it does not
 * use the actual vertex ids to determine the solution.
 *
 * Once a solution has been found the unique order of the vertices will be
 * returned
 */
class BranchSequenceLeaf {
 private:
  friend struct BranchSequenceLeafComparator;

  /// Why use a list, well because it has a lower memory overhead because
  /// leafs in the list can be removed when not in use, this is unlike a
  /// vector which is contiguous in memory,
  std::list<Index> getCanonicalVertexSequence_() const;

  bool canOrder_(const BranchSequenceLeaf& leaf) const;
  bool canOrder_(std::shared_ptr<BranchSequenceLeaf>& leaf) const;

  /**
   * @brief Checks to ensure that a branch is unique to a leaf
   *
   * A branch can be added to more than one leaf, but a single leaf
   * cannot contain the same branch more than once. It also cannot
   * add its self
   *
   * @return
   */
  bool branchIsUniqueToLeaf_(const Index branch_id) const;

  // sort by looking at the branches
  void sortLeaves_();
  void buildLabel_();

  /// Check if the branch exists in the current branch sequence
  // bool branchExists_(const Branch& branch) const;

  // Used to ensure branch is not looping around on itself
  Level level_;

  Index parallel_order_num_ = 0;

  // Separate from the branch id this is used to cluster groups of branches
  // together
  Index branch_group_id_;

  // The id of the branch
  Index branch_id_;

  // The actual branch
  Branch branch_;

  ContentLabel label_;

  // Pointers to leafs connected to 'branch_'
  std::vector<std::shared_ptr<BranchSequenceLeaf>> branch_sequence_;

  // Indicates whether the current branch leaf is ordered correctly
  // in the parent list
  bool sorted_ = false;

  // When a branch is the terminating branch this flag is set to true
  bool branch_end_ = false;

 public:
  BranchSequenceLeaf(const Index& branch_id,
                     const std::vector<Index>& branch_vertices,
                     const Index starting_vertex, const Graph& graph)
      : branch_id_(branch_id),
        branch_(branch_vertices, starting_vertex, graph) {
    std::cout << "Error 1" << std::endl;
    buildLabel_();
  };

  BranchSequenceLeaf(const Index& branch_id, const Branch& branch)
      : branch_id_(branch_id), branch_(branch) {
    std::cout << "Error 2" << std::endl;
    buildLabel_();
  };

  // Deep copy
  BranchSequenceLeaf(const BranchSequenceLeaf& leaf)
      : branch_id_(leaf.branch_id_),
        branch_(leaf.branch_),
        label_(leaf.label_),
        level_(leaf.level_),
        sorted_(leaf.sorted_),
        branch_end_(leaf.branch_end_) {
    for (const std::shared_ptr<BranchSequenceLeaf>& br_leaf :
         leaf.branch_sequence_) {
      branch_sequence_.push_back(std::shared_ptr<BranchSequenceLeaf>(br_leaf));
    }
    buildLabel_();
  };

  // Move constructor
  BranchSequenceLeaf(BranchSequenceLeaf&& leaf)
      : branch_id_(leaf.branch_id_),
        branch_(leaf.branch_),
        label_(leaf.label_),
        level_(leaf.level_),
        sorted_(leaf.sorted_),
        branch_end_(leaf.branch_end_) {

    branch_sequence_.insert(
        branch_sequence_.end(),
        std::make_move_iterator(leaf.branch_sequence_.begin()),
        std::make_move_iterator(leaf.branch_sequence_.end()));

    leaf.branch_sequence_.erase(leaf.branch_sequence_.begin(),
                                leaf.branch_sequence_.end());
  };

  Level getBranchLevel() const { return level_; }

  votca::Index getId() const { return branch_id_; }
  /// Does what it says, to do this it must parse the whole tree
  std::vector<Index> getCanonicalVertexSequence() const;

  /**
   * @brief build a branch string id that is based on the cummulitive contents
   * of all the branches in the sequence
   *
   * @return
   */
  ContentLabel getTreeContentLabel() const;

  bool isSorted() const noexcept { return sorted_; }

  /**
   * @brief Each branch sequence leaf should contain a sequence of at least
   * two branches or be a terminating branch, if neither condition is true
   * returns true (for incomplete)
   *
   * The premise behind this method is that each branch source and terminal
   * should be at a junction unless it is a terminating branch.
   *
   * @param
   *
   * @return
   */
  bool sequenceIsIncomplete() const noexcept {
    if (isDangling() || branch_sequence_.size() == 1) return true;
    return false;
  }

  /**
   * @brief Check to see if the tree as a whole is incomplete
   *
   * @return
   */
  bool treeIsIncomplete() const noexcept {
    if (sequenceIsIncomplete()) return true;
    for (const std::shared_ptr<BranchSequenceLeaf>& leaf_ : branch_sequence_) {
      if (leaf_->treeIsIncomplete()) return true;
    }
    return false;
  }

  void makeBranchEnd() {
    assert(branch_sequence_.size() == 0 &&
           "Cannot make branch end, there are branches in the sequence");
    branch_end_ = true;
  }

  /**
   * @brief If the branch sequence leaf has not been made a terminating leaf
   * but no branches have been added to the sequence will return true
   *
   * @return
   */
  bool isDangling() const noexcept {
    if (branch_end_ == false && branch_sequence_.size() == 0) return true;
    return false;
  }

  /**
   * @brief Check if the branches can be ordered by looking to see if
   * there exist any dangling branches
   *
   * As in branches that have not been specified as terminating but contain
   * no branches.
   *
   * @return
   */
  bool canOrder() const { return canOrder_(*this); }

  // Check if this leaf is the parent of the "leaf"
  // for this to be true this leaf's terminal vertex of its branch must
  // be the same as the source of "leaf"
  bool isRoot(std::shared_ptr<const BranchSequenceLeaf> leaf) const;

  // Order by looking at the branches
  void sortBranchSequence();

  /**
   * @brief Checks to see if the branch can be added to the sequence
   *
   * The source of 'branch' must match the terminal
   * of the internal 'branch_' and it must not have been previously added
   *
   * @param branch
   *
   * @return
   */
  // bool canAddBranch( const Index & branch_id, const Branch & branch) const;
  // Should be
  bool canAddLeaf(std::shared_ptr<const BranchSequenceLeaf> leaf) const;

  /**
   * @brief Adds a branch to the current branch sequence leaf the branch
   * can only be added if the "source" leaf of 'branch' is equal to
   * the end leaf of branch_.
   *
   * @param Branch
   */

  void addLeaf(std::shared_ptr<BranchSequenceLeaf> leaf);
};

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_BRANCHSEQUENCELEAF_H
