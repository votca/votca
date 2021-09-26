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

#ifndef VOTCA_CSG_PAIRLIST_H
#define VOTCA_CSG_PAIRLIST_H

// Standard includes
#include <map>
#include <vector>

// VOTCA includes
#include <votca/tools/types.h>

namespace votca {
namespace csg {

template <typename element_type, typename pair_type>
class PairList {
 public:
  PairList() = default;
  virtual ~PairList() { Cleanup(); }

  // this method takes ownership of p
  void AddPair(pair_type *p);

  using iterator = typename std::vector<pair_type *>::iterator;
  using const_iterator = typename std::vector<pair_type *>::const_iterator;
  typedef typename std::map<element_type, pair_type *> partners;

  iterator begin() { return pairs_.begin(); }
  iterator end() { return pairs_.end(); }

  const_iterator begin() const { return pairs_.begin(); }
  const_iterator end() const { return pairs_.end(); }
  pair_type *front() { return pairs_.front(); }
  pair_type *back() { return pairs_.back(); }
  bool empty() const { return pairs_.empty(); }

  Index size() const { return Index(pairs_.size()); }

  void Cleanup();

  pair_type *FindPair(element_type e1, element_type e2);

  const pair_type *FindPair(element_type e1, element_type e2) const;

  partners *FindPartners(element_type e1);

  using element_t = element_type;
  using pair_t = pair_type;

 protected:
  std::vector<pair_type *> pairs_;

  std::map<element_type, std::map<element_type, pair_type *>> pair_map_;
};

// this method takes ownership of p
template <typename element_type, typename pair_type>
inline void PairList<element_type, pair_type>::AddPair(pair_type *p) {
  /// \todo be careful, same pair object is used, some values might change (e.g.
  /// sign of distance vector)
  pair_map_[p->first()][p->second()] = p;
  pair_map_[p->second()][p->first()] = p;
  /// \todo check if unique
  pairs_.push_back(p);
}

template <typename element_type, typename pair_type>
inline void PairList<element_type, pair_type>::Cleanup() {
  for (auto &pair : pairs_) {
    delete pair;
  }
  pairs_.clear();
  pair_map_.clear();
}

template <typename element_type, typename pair_type>
inline pair_type *PairList<element_type, pair_type>::FindPair(element_type e1,
                                                              element_type e2) {
  typename std::map<element_type, std::map<element_type, pair_type *>>::iterator
      iter1;
  iter1 = pair_map_.find(e1);
  if (iter1 == pair_map_.end()) {
    return nullptr;
  }

  typename partners::iterator iter2;
  iter2 = iter1->second.find(e2);
  if (iter2 == iter1->second.end()) {
    return nullptr;
  }

  return iter2->second;
}

template <typename element_type, typename pair_type>
inline const pair_type *PairList<element_type, pair_type>::FindPair(
    element_type e1, element_type e2) const {
  typename std::map<element_type,
                    std::map<element_type, pair_type *>>::const_iterator iter1;
  iter1 = pair_map_.find(e1);
  if (iter1 == pair_map_.end()) {
    return nullptr;
  }

  typename partners::const_iterator iter2;
  iter2 = iter1->second.find(e2);
  if (iter2 == iter1->second.end()) {
    return nullptr;
  }

  return iter2->second;
}

template <typename element_type, typename pair_type>
typename PairList<element_type, pair_type>::partners *
    PairList<element_type, pair_type>::FindPartners(element_type e1) {
  typename std::map<element_type, std::map<element_type, pair_type *>>::iterator
      iter;
  if ((iter = pair_map_.find(e1)) == pair_map_.end()) {
    return nullptr;
  }
  return &(iter->second);
}

}  // namespace csg
}  // namespace votca

#endif /*  VOTCA_CSG_PAIRLIST_H */
