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

#include "branchtree.h"
#include "votca/tools/contentlabel.h"

namespace votca {
namespace tools {

void BranchTree::stepLevel() {
  assert(current_level_.num != 0 && "Cannot step level at source");
  current_level_.num--;
  prev_leafs_ = std::move(current_leafs_);
  current_leafs_.clear();
}

size_t BranchTree::countLeaves(Level level) const {
  if (level_ptrs_.count(level) == 0) return 0;
  return level_ptrs_.at(level).size();
}

void BranchTree::addBranch(votca::Index branch_id, Branch branch) {
  // Only check the source vertex when looking to add nodes
  assert(branch.getLevel() == current_level_ &&
         "Cannot add branch, wrong level");

  assert(branch_ids_.count(branch_id) == 0 &&
         "Cannot add branch it has already been added");
  auto new_leaf = std::shared_ptr<BranchSequenceLeaf>(
      new BranchSequenceLeaf(branch_id, branch));

  // Check to see if node is attached connected to any of the previously
  // added BranchSequenceNodes
  bool is_branch_end = true;
  for (std::shared_ptr<BranchSequenceLeaf>& leaf : prev_leafs_) {
    if (new_leaf->isRoot(leaf)) {
      new_leaf->addLeaf(leaf);
      is_branch_end = false;
    }
  }

  if (is_branch_end) new_leaf->makeBranchEnd();
  level_ptrs_[branch.getLevel()].push_back(new_leaf);
  current_leafs_.push_back(new_leaf);
}

std::vector<votca::Index> BranchTree::getBranchIds(Level level) const {
  std::vector<votca::Index> branch_ids;
  for (auto& leaf : level_ptrs_.at(level)) {
    branch_ids.push_back(leaf->getId());
  }
  return branch_ids;
}

// At each level we need to be able to calculate the string id
ContentLabel BranchTree::getContentLabel(Level level, votca::Index branch_id) {
  assert(level.num > -1 && "Cannot access branch in negative levels");
  assert(level >= current_level_ &&
         "Cannot access level, level has not been added to the tree");

  // Find the correct branch
  for (auto& leaf : level_ptrs_[level]) {
    if (leaf->getId() == branch_id) {
      if (not leaf->isSorted()) {
        leaf->sortBranchSequence();
      }
      return leaf->getTreeContentLabel();
    }
  }
  return ContentLabel();
}

}  // namespace tools
}  // namespace votca
