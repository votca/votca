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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Standard includes
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

// Third party includes
#include <boost/lexical_cast.hpp>

// Local VOTCA includes
#include "votca/tools/contentlabel.h"
#include "votca/tools/types.h"

using namespace boost;

namespace votca {
namespace tools {



/*
static vector<string> buildKeys_(const unordered_map<string, boost::any> vals) {
  vector<string> keys_temp;
  for (auto it : vals) {
    keys_temp.push_back(it.first);
  }
  sort(keys_temp.begin(), keys_temp.end());
  vector<std::string> keys;
  for (auto key_t : keys_temp) {
    auto it = vals.find(key);
    std::string key;
    key.append(key_t);
    key.append("=");
    keys.push_back(key);
  }
  return keys;
}

static vector<string> buildValues_(
    const unordered_map<string, boost::any> vals) {
  vector<string> keys;
  for (auto it : vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  vector<string> labels;
  for (auto key : keys) {
    auto it = vals.find(key);
    std::string label = "";
    label.append(convertToString_(it));
    label.append(",");
    labels.push_back(label);
  }
  return labels;
}*/

bool ContentLabel::isBranch_() const { //std::list<std::vector<KeyValType>> labels) const {
  if( type_ == BranchType::branch) return true;
//  if (labels.size() == 0) return false;
//  if (labels.front().size() == 0) return false;
//  if (labels.front().front()[0] == "{{") return true;
  return false;
}

/*
void ContentLabel::add(GraphNode gn) {
  if (labels_.back()[3] == "}") {
    throw std::runtime_argument(
        "Can only add graph node to content label before it has been made into "
        "a branch")
  }
  ContentLabel label = gn.getContentLabel();

  full_label_.insert(full_label_.end(), label.full_label_.begin(),
                     label.full_label_.end());
  brief_label_.insert(brief_label_.end(), label.brief_label_.begin(),
                      label.brief_label_.end());
}
void ContentLabel::add(Branch br) {

  if (full_label_.back().back() != "}") {
    throw std::runtime_argument(
        "Can only add branches to content label after it has been made into a "
        "branch")
  }
  full_label_.insert(full_label_.end(), label.full_label_.begin(),
                     label.full_label_.end());
  brief_label_.insert(brief_label_.end(), label.brief_label_.begin(),
                      label.brief_label_.end());

}*/

void ContentLabel::append(ContentLabel label) {

  // Check if current branch is stem and if label is furcation
  // First check if both content labels are considered branch labels
  bool this_label_is_branch = isBranch();
  bool that_label_is_branch = label.isBranch();

  if( type_ == BranchType::branch_stem && 
      label.type_ == BranchType::branch_furcation){
    std::cout << "******************" << std::endl;
    std::cout << "Appending branch furcation to stem" << std::endl;
    label_char_len_ += label.label_char_len_;
    BaseContentLabel::append_(label.labels_); 
    type_ = BranchType::tree;
  } else if( type_ == BranchType::tree &&
      label.type_ == BranchType::branch_furcation){
    std::cout << "******************" << std::endl;
    std::cout << "Appending branch furcation to tree" << std::endl;

    // On stem
    std::cout << "Back of tree should be ]" << std::endl;
    std::cout << labels_.back().back()[2] << labels_.back().back()[3] << std::endl;
    labels_.back().back()[3] = "),"; // remove ending ] and replace with ,
    std::cout << "After switch should be )," << std::endl;
    std::cout << labels_.back().back()[2] << labels_.back().back()[3] << std::endl;
    //  
    std::cout << "front of furcation should be [{ " << std::endl;
    std::cout << label.labels_.front().front()[0] << label.labels_.front().front()[1] << std::endl;
    label.labels_.front().front()[0] = "{"; // remove leading [  
    std::cout << "After switch should be {" << std::endl;
    std::cout << label.labels_.front().front()[0] << label.labels_.front().front()[1] << std::endl;
    label_char_len_ += label.label_char_len_;
    --label_char_len_;
    BaseContentLabel::append_(label.labels_); 

  } else if( this_label_is_branch && that_label_is_branch ) {
    std::runtime_error("Appending branch labels is not yet supported.");
  } else if (not this_label_is_branch && not that_label_is_branch) {
    // If it is not just append as normal sequentially
    label_char_len_ += label.label_char_len_;
    BaseContentLabel::append_(label.labels_); 
  } else {
    // This is an error cannot append content labels that are not both either
    // branch labels or are both not branch labels
    std::runtime_error("Cannot append, combination is not allowed.");
  }
}

bool ContentLabel::isBranch() const {
  return isBranch_();
} 

void ContentLabel::makeFurcationBranch() {
  if (not isBranch() ){
    std::string error_msg = "Cannot make into furcation stem must first be ";
    error_msg += "considered a branch.";
    throw std::runtime_error(error_msg);
  }

  auto front_node = labels_.front();
 
  size_t chars = 0; 
  for( auto values : front_node ) {
    for ( auto value : values ) {
      std::cout << "Value in front node " << value << std::endl;
      chars += value.size();
    } 
  }

  labels_.pop_front();
  if( labels_.size() == 0){
    return;
  }
  if( labels_.size() == 1){
    labels_.front().front()[0] = "[(";
    labels_.back().back()[3] = ")]";
  }else if(labels_.size() == 2) {
    labels_.front().front()[0] = "[{";
    labels_.back().back()[3] = ")]";
  }

  type_ = BranchType::branch_furcation;
  label_char_len_ += 2 - chars;
}

void ContentLabel::makeFurcationStem() {
  if (not isBranch() ){
    std::string error_msg = "Cannot make into furcation branch must first be ";
    error_msg += "considered a branch.";
    throw std::runtime_error(error_msg);
  }

  labels_.back().back()[3] = ")<";

  type_ = BranchType::branch_stem;

  ++label_char_len_;
}


void ContentLabel::makeBranchLabel() {
  std::cout << "calling contains branch" << std::endl;
  if (isBranch()) return;
  std::cout << "calling labels size" << std::endl;
  if (labels_.size() == 0) {
    std::vector<KeyValType> val;
    val.resize(1);
    labels_.push_back(val);
  }
  std::cout << "Assigning branch brackets" << std::endl;

  if( labels_.size() == 1 ) {
    if( labels_.front().size() == 0 ) {
      KeyValType val;
      labels_.front().push_back(val);
    }
    labels_.front().front()[0] = "(";
    // One of the } braces will replace a ;
    labels_.back().back()[3] = ")";
    label_char_len_ += 1;
  } else {

    if( labels_.front().size() == 0 ) {
      KeyValType val;
      labels_.front().push_back(val);
    }
    if( labels_.back().size() == 0 ) {
      KeyValType val;
      labels_.back().push_back(val);
    }
    labels_.front().front()[0] = "(";
    labels_.front().back()[3] = ")";
    labels_.back().front()[0] = "(";
    labels_.back().back()[3] = ")";
    // Two of the } will replace a semicolon ;
    label_char_len_ += 2;
  
    if ( labels_.size() > 2) {
      std::list<std::vector<KeyValType>>::iterator iter = labels_.end();
      --iter;
      --iter;
      // Second to last element remove the ; 
      iter->back()[3] = "}";

      iter = labels_.begin();
      ++iter;
      iter->front()[0] = "{";
      ++label_char_len_;
    }

  }
  type_ = BranchType::branch;
}

void ContentLabel::reverse() {
  if( isBranch() ) {
    // The order is important do not change
    auto iter = labels_.end();
    --iter;
    iter->back()[3] = "}"; 
    if ( labels_.size() > 3) {
      --iter;
      iter->back()[3] = ";"; 
    }
    iter = labels_.begin();
    iter->front()[0] = "{";
    iter->back()[3] = "}}";
    if ( labels_.size() > 3) {
      ++iter;
      iter->back()[3] = ""; 
    }
    iter = labels_.end();
    --iter;
    iter->front()[0] = "{{"; 
  }
  BaseContentLabel::reverse();
}

void ContentLabel::calcCharLen_() {
  label_char_len_ = get().length();
}

}  // namespace tools
}  // namespace votca
