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

namespace votca {
namespace tools {
/*
class BranchTree {
  private:
    std::map<Index, std::vector<BranchSequenceNode *>> level_ptrs;

    int current_level_;

    std::list<BranchSequenceNode *> prev_nodes
    std::list<BranchSequenceNode *> current_nodes;
  public:
    BranchTree(int levels) : level_ptrs(levels), current_level_(levels) {};

    void stepLevel(){
      assert(level!=0 && "Cannot step level at source");

      level--;
      prev_nodes = std::move(current_nodes);
      current_nodes.clear();
    }

    void addBranch(int level, Branch branch) {
      // Only check the source vertex when loolking to add nodes
      assert(level==current_level_ && "Cannot add branch, wrong level");


      BranchSequenceNode* node = new BranchSequenceNode(branch);

      // Check to see if node is attached connected to any of the previously
      // added BranchSequenceNodes
      for ( BranchSequenceNode * p_node : prev_nodes) {
        if( node->isParent(*p_node) ){
          node->addLeaf(p_node);
        }
      }
      current_nodes.push_back(node)

    }

    // At each level we need to be able to calculate the string id
    std::string getStrId(int level, int branch_id ){
      assert(level>=0 && "Cannot access branch in negative levels");
      assert(level>=current_level_ && "Cannot access level, level has not been
added to the tree");
      // Find the correct branch
      for ( BranchSequenceNode * node : level_ptrs[level]){
        if(node->branch_id_ == branch_id) {
          return node->getCummulitiveBranchStringId();
        }
      }
      return "";
    }
};
*/
}
}  // namespace votca
#endif  // VOTCA_TOOLS_BRANCHTREE_H
