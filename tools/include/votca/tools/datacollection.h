/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_TOOLS_DATACOLLECTION_H
#define VOTCA_TOOLS_DATACOLLECTION_H

// Standard includes
#include <cassert>
#include <map>
#include <sstream>
#include <vector>

// Local VOTCA includes
#include "tokenizer.h"
#include "types.h"

namespace votca {
namespace tools {

/**
 * \brief This class handles a set of arrays which can be identified by name
 * tags
 *
 * This class is a Container of arrays. The arrays can be accessed by
 * specifying a name or wildcard.
 *
 **/
template <typename T>
class DataCollection {
 public:
  /**
   * \brief The array class, extends vector by a name tag
   */
  class array : public std::vector<T> {
   public:
    array(std::string name) { name_ = name; }
    const std::string &getName() const { return name_; }

   private:
    std::string name_;
  };

  using iterator = typename std::vector<array *>::iterator;
  using const_iterator = typename std::vector<array *>::const_iterator;

  /**
   * \brief class for array selection
   */
  class selection {
   public:
    selection() = default;
    ~selection() = default;

    Index size() const { return Index(arrays_.size()); }
    bool empty() const { return arrays_.empty(); }
    array &operator[](Index i) {
      assert(i < Index(arrays_.size()));
      return *(arrays_[i]);
    }

    const array &operator[](Index i) const {
      assert(i < Index(arrays_.size()));
      return *(arrays_[i]);
    }

    void push_back(array *a) { arrays_.push_back(a); }
    void push_back(selection *s) {
      arrays_.insert(arrays_.end(), s->begin(), s->end());
    }

    iterator begin() { return arrays_.begin(); }
    iterator end() { return arrays_.end(); }
    const_iterator begin() const { return arrays_.begin(); }
    const_iterator end() const { return arrays_.end(); }

   private:
    std::vector<array *> arrays_;
  };

  /// constructor
  DataCollection() = default;
  /// destructor
  ~DataCollection() { clear(); }

  /**
   * \brief clears the data collection
   */
  void clear();
  /**
   *  \ brief returns the number of arrays
   */
  Index size() const { return Index(data_.size()); }
  bool empty() const { return data_.empty(); }
  array &operator[](Index i) {
    assert(i < Index(data_.size()));
    return *(data_[i]);
  }
  iterator begin() { return data_.begin(); }
  iterator end() { return data_.end(); }
  const_iterator begin() const { return data_.begin(); }
  const_iterator end() const { return data_.end(); }

  /**
   * \brief create a new array
   */
  array *CreateArray(std::string name);

  /**
   * \brief access the data container
   */
  std::vector<array *> &Data() { return data_; }
  const std::vector<array *> &Data() const { return data_; }

  /**
   * \brief access an array by name
   */
  array *ArrayByName(std::string name);

  /**
   * \brief select a set of arrays
   *
   * WARNING If attempting to append to an existing selection you must be
   * careful if there exist more than one array with the same name the
   * first array name that matches 'strselection' will be appended.
   */
  selection *select(std::string strselection, selection *sel_append = nullptr);

 private:
  std::vector<array *> data_;

  std::map<std::string, array *> array_by_name_;
};

template <typename T>
void DataCollection<T>::clear() {

  for (auto &d : data_) {
    delete d;
  }
  data_.clear();
}

template <typename T>
typename DataCollection<T>::array *DataCollection<T>::CreateArray(
    std::string name) {
  assert(ArrayByName(name) == nullptr);
  array *a = new array(name);
  data_.push_back(a);
  array_by_name_[name] = a;

  return a;
}

template <typename T>
typename DataCollection<T>::array *DataCollection<T>::ArrayByName(
    std::string name) {
  typename std::map<std::string, array *>::iterator i =
      array_by_name_.find(name);
  if (i == array_by_name_.end()) {
    return nullptr;
  }
  return (*i).second;
}

template <typename T>
typename DataCollection<T>::selection *DataCollection<T>::select(
    std::string strselection, selection *sel_append) {

  typename DataCollection<T>::selection *sel = sel_append;
  if (!sel_append) {
    sel = new typename DataCollection<T>::selection;
  }

  for (auto &pair : array_by_name_) {
    if (wildcmp(strselection.c_str(), pair.second->getName().c_str())) {
      sel->push_back(pair.second);
    }
  }
  return sel;
}

std::ostream &operator<<(std::ostream &out,
                         const DataCollection<double>::selection &sel);
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_DATACOLLECTION_H
