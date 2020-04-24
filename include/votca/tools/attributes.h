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

#ifndef VOTCA_TOOLS_ATTRIBUTES_H
#define VOTCA_TOOLS_ATTRIBUTES_H
#pragma once

#include "votca/tools/types.h"
#include <boost/any.hpp>
#include <string>
#include <unordered_map>

namespace votca {
namespace tools {

enum class LabelType { verbose, concise };

class Attributes {
 private:
  std::string full_label_{""};
  std::string brief_label_{""};

  // Consists of a key and a value
  std::unordered_map<std::string, boost::any> attributes_;
  void buildLabels_();

 public:
  Attributes() = default;

  /// Constructor
  Attributes(const std::unordered_map<std::string, boost::any> values);

  void set(const std::unordered_map<std::string, boost::any> values);
  void reset(const std::string key, const boost::any value);
  void add(const std::string key, const boost::any value);
  bool exists(const std::string& key) const;

  template <typename T>
  T get(const std::string key) const {
    if (attributes_.count(key) == 0)
      throw std::invalid_argument("Cannot get attribute key does not exist");
    T val;
    try {
      val = boost::any_cast<T>(attributes_.at(key));
    } catch (const boost::bad_any_cast& e) {
      throw boost::bad_any_cast("Bad any cast in attributes class");
    }
    return val;
  }

  std::string getContentLabel(LabelType type = LabelType::verbose) const {
    if (type == LabelType::verbose) return full_label_;
    return brief_label_;
  }
  bool operator==(const Attributes attr) const;
  bool operator!=(const Attributes attr) const;

  friend std::ostream& operator<<(std::ostream& os, const Attributes attr);
};

bool cmpAttributes(Attributes attr1, Attributes attr2);

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_ATTRIBUTES_H
