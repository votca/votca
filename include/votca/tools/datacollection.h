/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_TOOLS_DATACOLLECTION_H
#define _VOTCA_TOOLS_DATACOLLECTION_H

#include "tokenizer.h"
#include "types.h"
#include <assert.h>
#include <map>
#include <sstream>
#include <vector>

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
    array(std::string name) { _name = name; }
    const std::string &getName() const { return _name; }

   private:
    std::string _name;
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

    Index size() const { return Index(_arrays.size()); }
    bool empty() const { return _arrays.empty(); }
    array &operator[](Index i) {
      assert(i < Index(_arrays.size()));
      return *(_arrays[i]);
    }

    const array &operator[](Index i) const {
      assert(i < Index(_arrays.size()));
      return *(_arrays[i]);
    }

    void push_back(array *a) { _arrays.push_back(a); }
    void push_back(selection *s) {
      _arrays.insert(_arrays.end(), s->begin(), s->end());
    }

    iterator begin() { return _arrays.begin(); }
    iterator end() { return _arrays.end(); }
    const_iterator begin() const { return _arrays.begin(); }
    const_iterator end() const { return _arrays.end(); }

   private:
    std::vector<array *> _arrays;
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
  Index size() const { return Index(_data.size()); }
  bool empty() const { return _data.empty(); }
  array &operator[](Index i) {
    assert(i < Index(_data.size()));
    return *(_data[i]);
  }
  iterator begin() { return _data.begin(); }
  iterator end() { return _data.end(); }
  const_iterator begin() const { return _data.begin(); }
  const_iterator end() const { return _data.end(); }

  /**
   * \brief create a new array
   */
  array *CreateArray(std::string name);

  /**
   * \brief access the data container
   */
  std::vector<array *> &Data() { return _data; }
  const std::vector<array *> &Data() const { return _data; }

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
  std::vector<array *> _data;

  std::map<std::string, array *> _array_by_name;
};

template <typename T>
void DataCollection<T>::clear() {

  for (auto &d : _data) {
    delete d;
  }
  _data.clear();
}

template <typename T>
typename DataCollection<T>::array *DataCollection<T>::CreateArray(
    std::string name) {
  assert(ArrayByName(name) == nullptr);
  array *a = new array(name);
  _data.push_back(a);
  _array_by_name[name] = a;

  return a;
}

template <typename T>
typename DataCollection<T>::array *DataCollection<T>::ArrayByName(
    std::string name) {
  typename std::map<std::string, array *>::iterator i =
      _array_by_name.find(name);
  if (i == _array_by_name.end()) {
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

  for (auto &pair : _array_by_name) {
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

#endif  // _VOTCA_TOOLS_DATACOLLECTION_H
