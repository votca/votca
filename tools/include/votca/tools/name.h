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

#ifndef VOTCA_TOOLS_NAME_H
#define VOTCA_TOOLS_NAME_H

// Standard includes
#include <cassert>
#include <string>

namespace votca {
namespace tools {

/**
 * \brief Name object
 *
 * This object is meant to be used a derived type for larger objects. In this
 * way the same methods for setting and getting the name of an object will be
 * uniformly defined.
 *
 */
class Name {
 private:
  std::string name_;
  bool name_set_{false};

 public:
  Name() = default;
  Name(const std::string name) : name_(name), name_set_(true) {};
  void setName(const std::string &name) {
    name_ = name;
    name_set_ = true;
  }
  const std::string &getName() const {
    assert(name_set_ && "No name has been set, cannot get name.");
    return name_;
  }
};
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_NAME_H
