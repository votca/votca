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

#ifndef VOTCA_TOOLS_CONTENTLABEL_H
#define VOTCA_TOOLS_CONTENTLABEL_H
#pragma once

// Standard includes
#include <array>
#include <list>
#include <string>
#include <unordered_map>
#include <vector>

// Third party includes
#include <boost/any.hpp>

// Local VOTCA includes
#include "contentlabelbase.h"

namespace votca {
namespace tools {

class GraphNode;
class Branch;

class ContentLabel : public BaseContentLabel {
  private:
    /// Used to update the char length should avoid using if possible
    void calcCharLen_();

    // Determine based on labels if it is a branch
    bool isBranch_(std::list<std::vector<KeyValType>> labels) const;
 public:
  ContentLabel() = default;
  ContentLabel(std::unordered_map<std::string, boost::any> values) :  
    BaseContentLabel(values) {};

  // void add(GraphNode gn);
  //void add(Branch br);
  bool isBranch() const;
  
  /**
   * @brief if the content label has been made a branch label only content 
   * labels that are also branches can be appended
   *
   * @param ContentLabel
   */
  virtual void append(ContentLabel) final;
 
  /**
   * @brief When a content label is made a branch curly braces are appended
   * to either side of the content label, and an extra set of curly braces
   * are placed around the end nodes
   */
  void makeBranchLabel();

};

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_CONTENTLABEL_H
