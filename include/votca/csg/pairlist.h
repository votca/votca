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

#ifndef _VOTCA_CSG_PAIRLIST_H
#define _VOTCA_CSG_PAIRLIST_H

#include <map>
#include <vector>

namespace votca {
namespace csg {

template <typename element_type, typename pair_type>
class PairList {
 public:
  PairList() {}
  virtual ~PairList() { Cleanup(); }

  // this method takes ownership of p
  void AddPair(pair_type *p);

  typedef typename std::vector<pair_type *>::iterator iterator;
  typedef typename std::vector<pair_type *>::const_iterator const_iterator;
  typedef typename std::map<element_type, pair_type *> partners;

  iterator begin() { return _pairs.begin(); }
  iterator end() { return _pairs.end(); }

  const_iterator begin() const { return _pairs.begin(); }
  const_iterator end() const { return _pairs.end(); }
  pair_type *front() { return _pairs.front(); }
  pair_type *back() { return _pairs.back(); }
  bool empty() const { return _pairs.empty(); }

  int size() const { return _pairs.size(); }

  void Cleanup();

  pair_type *FindPair(element_type e1, element_type e2);

  const pair_type *FindPair(element_type e1, element_type e2) const;

  partners *FindPartners(element_type e1);

  typedef element_type element_t;
  typedef pair_type pair_t;

 protected:
  std::vector<pair_type *> _pairs;

  std::map<element_type, std::map<element_type, pair_type *>> _pair_map;
};

// this method takes ownership of p
template <typename element_type, typename pair_type>
inline void PairList<element_type, pair_type>::AddPair(pair_type *p) {
  /// \todo be careful, same pair object is used, some values might change (e.g.
  /// sign of distance vector)
  _pair_map[p->first()][p->second()] = p;
  _pair_map[p->second()][p->first()] = p;
  /// \todo check if unique
  _pairs.push_back(p);
}

template <typename element_type, typename pair_type>
inline void PairList<element_type, pair_type>::Cleanup() {
  for (iterator iter = _pairs.begin(); iter != _pairs.end(); ++iter)
    delete *iter;
  _pairs.clear();
  _pair_map.clear();
}

template <typename element_type, typename pair_type>
inline pair_type *PairList<element_type, pair_type>::FindPair(element_type e1,
                                                              element_type e2) {
  typename std::map<element_type, std::map<element_type, pair_type *>>::iterator
      iter1;
  iter1 = _pair_map.find(e1);
  if (iter1 == _pair_map.end()) return nullptr;

  typename partners::iterator iter2;
  iter2 = iter1->second.find(e2);
  if (iter2 == iter1->second.end()) return nullptr;

  return iter2->second;
}

template <typename element_type, typename pair_type>
inline const pair_type *PairList<element_type, pair_type>::FindPair(
    element_type e1, element_type e2) const {
  typename std::map<element_type,
                    std::map<element_type, pair_type *>>::const_iterator iter1;
  iter1 = _pair_map.find(e1);
  if (iter1 == _pair_map.end()) return nullptr;

  typename partners::const_iterator iter2;
  iter2 = iter1->second.find(e2);
  if (iter2 == iter1->second.end()) return nullptr;

  return iter2->second;
}

template <typename element_type, typename pair_type>
typename PairList<element_type, pair_type>::partners *
    PairList<element_type, pair_type>::FindPartners(element_type e1) {
  typename std::map<element_type, std::map<element_type, pair_type *>>::iterator
      iter;
  if ((iter = _pair_map.find(e1)) == _pair_map.end()) return nullptr;
  return &(iter->second);
}

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_PAIRLIST_H */
