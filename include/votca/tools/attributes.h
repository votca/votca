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
  void checkKey_(const std::string& key);

 public:
  Attributes() = default;

  /// Constructor
  //  Attributes(const std::unordered_map<std::string, boost::any> values);

  template <typename T>
  Attributes(const std::unordered_map<std::string, T> values) {
    for (auto key_val : values) attributes_[key_val.first] = key_val.second;
    buildLabels_();
  }

  template <typename T, typename U>
  Attributes(const std::unordered_map<std::string, T> values1,
             const std::unordered_map<std::string, U> values2) {
    for (auto key_val : values1) attributes_[key_val.first] = key_val.second;
    for (auto key_val : values2) attributes_[key_val.first] = key_val.second;
    buildLabels_();
  }

  template <typename T, typename U, typename V>
  Attributes(const std::unordered_map<std::string, T> values1,
             const std::unordered_map<std::string, U> values2,
             const std::unordered_map<std::string, V> values3) {
    for (auto key_val : values1) attributes_[key_val.first] = key_val.second;
    for (auto key_val : values2) attributes_[key_val.first] = key_val.second;
    for (auto key_val : values3) attributes_[key_val.first] = key_val.second;
    buildLabels_();
  }
  /// Returns a list of the keys
  std::vector<std::string> getKeys() const noexcept;

  /// Returns the keys and the of the values that they store
  std::unordered_map<std::string, std::string> getTypes() const noexcept;
  /// Sets the attributes with the provided values
  //  void set(const std::unordered_map<std::string, boost::any> values);

  template <typename T>
  void set(const std::unordered_map<std::string, T> values) {
    for (auto& key_val : values) {
      checkKey_(key_val.first);
    }
    attributes_.clear();
    for (auto& key_val : values) {
      attributes_[key_val.first] = key_val.second;
    }
    buildLabels_();
  }

  /// Resets the attributes with the provided value, the key must already exist
  template <typename T>
  void reset(const std::string key, const T value) {
    checkKey_(key);
    if (attributes_.count(key) == 0) {
      throw std::invalid_argument("Cannot reset attribute, " + key +
                                  " the key is not known.");
    }
    if (attributes_.at(key).type() != typeid(T)) {
      throw std::invalid_argument(
          "Cannot reset attribute to a different type.");
    }
    attributes_[key] = value;
    buildLabels_();
  }

  template <typename T>
  void add(const std::string key, const T value) {
    assert(attributes_.count(key) == 0 &&
           "Cannot add attribute Attributes instance, the key has already been "
           "used");
    checkKey_(key);
    attributes_[key] = value;
    buildLabels_();
  }

  template <typename T>
  void add(const std::unordered_map<std::string, T> values) {
    for (auto& key_val : values) {
      checkKey_(key_val.first);
      assert(
          attributes_.count(key_val.first) == 0 &&
          "Cannot add attribute Attributes instance, the key has already been "
          "used");
      attributes_[key_val.first] = key_val.second;
    }
    buildLabels_();
  }

  void remove(const std::string& key);

  template <typename T>
  void removeAll() {
    auto it = attributes_.begin();
    while (it != attributes_.end()) {
      if (it->second.type() == typeid(T)) {
        it = attributes_.erase(it);
      } else {
        ++it;
      }
    }
  }

  bool exists(const std::string& key) const;

  template <typename T>
  T get(const std::string key) const {
    if (attributes_.count(key) == 0)
      throw std::invalid_argument("Cannot get attribute key does not exist");
    T val;
    try {
      val = boost::any_cast<T>(attributes_.at(key));
    } catch (const boost::bad_any_cast& e) {
      throw std::runtime_error("Bad any cast in attributes class");
    }
    return val;
  }

  /// Return all values of a specific type
  template <typename T>
  std::unordered_map<std::string, T> getAll() const {
    std::unordered_map<std::string, T> values;
    for (auto& key_val : attributes_) {
      if (key_val.second.type() == typeid(T)) {
        values[key_val.first] = boost::any_cast<T>(key_val.second);
      }
    }
    return values;
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
