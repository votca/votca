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

#include "branchsequenceleaf.h"
#include "branch.h"
#include "votca/tools/contentlabel.h"
#include "votca/tools/graphnode.h"
#include "votca/tools/reducededge.h"
#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace votca {
namespace tools {

/**************************************************************************
 * Private methods
 **************************************************************************/
// Should sort the branches such that the shortest branch always comes first
// If branches are equal than they are sorted alpha numerically
struct BranchSequenceLeafComparator {
  inline bool operator()(std::shared_ptr<BranchSequenceLeaf>& n1,
                         std::shared_ptr<BranchSequenceLeaf>& n2) {
    if (n1->sortAndGetBranchSequenceIds_().length() >
        n2->sortAndGetBranchSequenceIds_().length())
      return false;
    if (n1->sortAndGetBranchSequenceIds_().length() <
        n2->sortAndGetBranchSequenceIds_().length())
      return true;
    return n1->getTreeContentLabel().compare(n2->getTreeContentLabel()) < 0;
  }
};

bool BranchSequenceLeaf::canOrder_(const BranchSequenceLeaf& leaf) const {
  if (sequenceIsIncomplete()) return false;
  for (std::shared_ptr<BranchSequenceLeaf> leaf_ : leaf.branch_sequence_) {
    if (leaf.canOrder_(leaf_) == false) return false;
  }
  return true;
}

bool BranchSequenceLeaf::canOrder_(
    std::shared_ptr<BranchSequenceLeaf>& leaf) const {
  if (sequenceIsIncomplete()) return false;
  for (std::shared_ptr<BranchSequenceLeaf> leaf_ : leaf->branch_sequence_) {
    if (leaf->canOrder_(leaf_) == false) return false;
  }
  return true;
}

bool BranchSequenceLeaf::branchExists_(const Branch& branch) const {
  for (const std::shared_ptr<BranchSequenceLeaf>& leaf : branch_sequence_) {
    if (leaf->branch_.getContentLabel().compare(branch.getContentLabel())) {
      return true;
    }
  }
  return false;
}

std::string BranchSequenceLeaf::sortAndGetBranchSequenceIds_() {
  //  c_branch_str_id_ = "{";
  //  c_branch_str_id_.append(branch_.getContentLabel());
  //  c_branch_str_id_.back() = '}';
  std::sort(branch_sequence_.begin(), branch_sequence_.end(),
            BranchSequenceLeafComparator());

  for (std::shared_ptr<BranchSequenceLeaf>& leaf_ : branch_sequence_) {
    c_branch_str_id_ += leaf_->c_branch_str_id_;
  }
  sorted_ = true;
  return c_branch_str_id_;
}

std::list<Index> BranchSequenceLeaf::getCanonicalVertexSequence_() const {
  std::vector<Index> seq_vec = branch_.getSequence();
  std::list<Index> vertex_sequence(seq_vec.begin(), seq_vec.end());
  for (const std::shared_ptr<BranchSequenceLeaf>& leaf_ : branch_sequence_) {
    std::list<Index> seq = leaf_->getCanonicalVertexSequence_();
    vertex_sequence.splice(vertex_sequence.end(), seq);
  }
  return vertex_sequence;
}

bool BranchSequenceLeaf::branchIsUniqueToLeaf_(const Index branch_id) const {
  if (branch_id_ == branch_id) return false;
  for (const std::shared_ptr<BranchSequenceLeaf>& leaf_ : branch_sequence_) {
    if (leaf_->branch_id_ == branch_id) return false;
  }
  return true;
}

bool BranchSequenceLeaf::canAddLeaf(
    std::shared_ptr<const BranchSequenceLeaf> leaf) const {
  return branchIsUniqueToLeaf_(leaf->branch_id_);
}

bool BranchSequenceLeaf::isRoot(
    std::shared_ptr<const BranchSequenceLeaf> leaf) const {
  if (branch_.getTerminal() == leaf->branch_.getSource()) return true;
  return false;
}
// Order by looking at the branches
void BranchSequenceLeaf::sortBranchSequence() {
  c_branch_str_id_ = sortAndGetBranchSequenceIds_();
}

std::string BranchSequenceLeaf::getTreeContentLabel() const {
  if (not sorted_)
    throw std::runtime_error(
        "Cannot get cummulitive branch string id because the branch sequence "
        "leaf has not been sorted");
  return c_branch_str_id_;
}

std::vector<Index> BranchSequenceLeaf::getCanonicalVertexSequence() const {
  std::list<Index> vertex_sequence = getCanonicalVertexSequence_();
  std::vector<Index> unique_vertex_sequence;
  std::unordered_set<Index> unique_vertices;
  while (vertex_sequence.size() > 0) {
    if (unique_vertices.count(vertex_sequence.front()) == 0) {
      unique_vertices.insert(vertex_sequence.front());
      unique_vertex_sequence.push_back(vertex_sequence.front());
    }
    vertex_sequence.pop_front();
  }
  return unique_vertex_sequence;
}

/*  void BranchSequenceLeaf::addBranch( const Index & branch_id, Branch branch)
  {

    if( canAddBranch(branch_id,branch)==false ) {
      throw std::runtime_error("Cannot add branch to branch sequence, it is
  either not unique or does not share a vertex");
    }
    /// Sequence needs to be reversed if the terminal vertex of 'branch_'
    /// does not match the source branch of 'branch' but does match the
    /// terminal vertex of 'branch'
    if ( branch_.getTerminal()==branch.getTerminal() &&
        branch_.getTerminal()!=branch.getSource()){
      branch.reverseSequence();
    }

    auto leaf = std::shared_ptr<BranchSequenceLeaf>(new
  BranchSequenceLeaf(branch_id, branch));
    branch_sequence_.push_back(std::move(leaf));
  }*/

void BranchSequenceLeaf::addLeaf(std::shared_ptr<BranchSequenceLeaf> leaf) {

  /// Check that the branch has not been previously added
  if (canAddLeaf(leaf) == false) {
    throw std::runtime_error(
        "Cannot add branch sequence leaf to branch sequence, it is either not "
        "unique or does not share a vertex");
  }

  /// Sequence needs to be reversed if the therminal vertex of 'branch_'
  /// does not match the source branch of 'branch' but does match the
  /// terminal vertex of 'branch'
  if (branch_.getTerminal() == leaf->branch_.getTerminal() &&
      branch_.getTerminal() != leaf->branch_.getSource()) {
    leaf->branch_.reverseSequence();
  }

  //    std::shared_ptr<BranchSequenceLeaf> leaf_ptr =
  //    std::shared_ptr<BranchSequenceLeaf>(new BranchSequenceLeaf(leaf));
  branch_sequence_.push_back(leaf);
}

}  // namespace tools
}  // namespace votca
