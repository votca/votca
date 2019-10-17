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

#ifndef _VOTCA_CSG_TRIPLELIST_H
#define _VOTCA_CSG_TRIPLELIST_H

#include <map>
#include <vector>

namespace votca {
namespace csg {

template <typename element_type, typename triple_type>
class TripleList {
 public:
  TripleList() {}
  virtual ~TripleList() { Cleanup(); }

  void AddTriple(triple_type *t);

  typedef typename std::vector<triple_type *>::iterator iterator;

  iterator begin() { return _triples.begin(); }
  iterator end() { return _triples.end(); }
  typename std::vector<triple_type *>::size_type size() {
    return _triples.size();
  }
  triple_type *front() { return _triples.front(); }
  triple_type *back() { return _triples.back(); }
  bool empty() { return _triples.empty(); }

  void Cleanup();

  triple_type *FindTriple(element_type e1, element_type e2, element_type e3);

  typedef element_type element_t;
  typedef triple_type triple_t;

 private:
  std::vector<triple_type *> _triples;

  std::map<element_type,
           std::map<element_type, std::map<element_type, triple_type *>>>
      _triple_map;
};

template <typename element_type, typename triple_type>
inline void TripleList<element_type, triple_type>::AddTriple(triple_type *t) {
  //(*t)[i] gives access to ith element of tuple object (i=0,1,2).
  // only consider the permutations of elements (1,2) of the tuple object ->
  // tuple objects of the form (*,1,2) and (*,2,1) are considered to be the same
  _triple_map[std::get<0>(*t)][std::get<1>(*t)][std::get<2>(*t)] = t;
  _triple_map[std::get<0>(*t)][std::get<2>(*t)][std::get<1>(*t)] = t;
  /// \todo check if unique
  _triples.push_back(t);
}

template <typename element_type, typename triple_type>
inline void TripleList<element_type, triple_type>::Cleanup() {
  for (iterator iter = _triples.begin(); iter != _triples.end(); ++iter)
    delete *iter;
  _triples.clear();
  _triple_map.clear();
}

template <typename element_type, typename triple_type>
inline triple_type *TripleList<element_type, triple_type>::FindTriple(
    element_type e1, element_type e2, element_type e3) {
  typename std::map<
      element_type,
      std::map<element_type, std::map<element_type, triple_type *>>>::iterator
      iter1;
  iter1 = _triple_map.find(e1);
  if (iter1 == _triple_map.end()) return nullptr;

  typename std::map<element_type,
                    std::map<element_type, triple_type *>>::iterator iter2;
  iter2 = iter1->second.find(e2);
  if (iter2 == iter1->second.end()) return nullptr;

  typename std::map<element_type, triple_type *>::iterator iter3;
  iter3 = iter2->second.find(e3);
  if (iter3 == iter2->second.end()) return nullptr;

  return iter3->second;
}

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_TRIPLELIST_H */
