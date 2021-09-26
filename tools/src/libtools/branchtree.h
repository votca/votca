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
#include "votca/tools/contentlabel.h"
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

  // each branch must have a unique id for this to work
  std::set<votca::Index> branch_ids_;
  std::list<std::shared_ptr<BranchSequenceLeaf>> prev_leafs_;
  std::list<std::shared_ptr<BranchSequenceLeaf>> current_leafs_;

 public:
  BranchTree() = default;
  BranchTree(Level max_level) : current_level_(max_level){};

  void stepLevel();

  Level getCurrentLevel() const { return current_level_; }

  size_t countLeaves(Level level) const;

  std::vector<votca::Index> getBranchIds(Level level) const;

  void addBranch(votca::Index branch_id, Branch branch);

  // At each level we need to be able to calculate the string id
  ContentLabel getContentLabel(Level level, votca::Index branch_id);
};

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_BRANCHTREE_H
