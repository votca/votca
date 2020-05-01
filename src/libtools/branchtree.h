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

#ifndef VOTCA_TOOLS_BRANCHTREE_H
#define VOTCA_TOOLS_BRANCHTREE_H
#pragma once

#include "branchsequenceleaf.h"
#include "votca/tools/types.h"

#include <map>
#include <memory>
#include <vector>

namespace votca {
namespace tools {

/**
 * @brief BranchTree
 *
 * It is responsible for all branch sequence leafs, memory for leafs can only
 * be created by the BranchTree
 */
class BranchTree {
 private:
  // The levels of the different branches
  typedef std::map<Level, std::vector<std::shared_ptr<BranchSequenceLeaf>>>
      LevelManager;
  LevelManager level_ptrs_;

  Level current_level_;

  std::list<std::shared_ptr<BranchSequenceLeaf>> prev_leafs_;
  std::list<std::shared_ptr<BranchSequenceLeaf>> current_leafs_;

 public:
  BranchTree() = default;
  BranchTree(Level max_level) : current_level_(max_level){};

  void stepLevel() {
    assert(current_level_.num != 0 && "Cannot step level at source");
    current_level_.num--;
    prev_leafs_ = std::move(current_leafs_);
    current_leafs_.clear();
  }

  Level getCurrentLevel() const { return current_level_; }

  size_t countLeaves(Level level) const {
    if (level_ptrs_.count(level) == 0) return 0;
    return level_ptrs_.at(level).size();
  }

  void addBranch(Level level, votca::Index branch_id, Branch branch) {
    // Only check the source vertex when looking to add nodes
    assert(level == current_level_ && "Cannot add branch, wrong level");

    auto new_leaf = std::shared_ptr<BranchSequenceLeaf>(
        new BranchSequenceLeaf(branch_id, branch));

    // Check to see if node is attached connected to any of the previously
    // added BranchSequenceNodes
    for (std::shared_ptr<BranchSequenceLeaf>& leaf : prev_leafs_) {
      if (new_leaf->isRoot(leaf)) {
        new_leaf->addLeaf(leaf);
      }
    }
    level_ptrs_[level].push_back(new_leaf);
    current_leafs_.push_back(new_leaf);
  }

  // At each level we need to be able to calculate the string id
  std::string getContentLabel(Level level, votca::Index branch_id) {
    assert(level > negative_level && "Cannot access branch in negative levels");
    assert(level >= current_level_ &&
           "Cannot access level, level has not been added to the tree");
    // Find the correct branch
    for (auto& leaf : level_ptrs_[level]) {
      if (leaf->getId() == branch_id) {
        return leaf->getTreeContentLabel();
      }
    }
    return "";
  }
};

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_BRANCHTREE_H
