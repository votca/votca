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

#ifndef VOTCA_CSG_EXCLUSIONLIST_H
#define VOTCA_CSG_EXCLUSIONLIST_H

// Standard includes
#include <iostream>
#include <list>
#include <map>

// Local VOTCA includes
#include "bead.h"

namespace votca {
namespace csg {

/// \todo fill _excl_by_bead
/// \todo no ids but pointers, use PairList

class Topology;
class Bead;

class ExclusionList {
 public:
  ExclusionList() = default;
  ~ExclusionList() { Clear(); }

  void Clear(void);

  template <typename iterable>
  void Remove(iterable &l);

  template <typename iterable>
  void ExcludeList(iterable &l);

  struct exclusion_t {
    Bead *_atom;
    std::list<Bead *> _exclude;
  };

  void CreateExclusions(Topology *top);
  exclusion_t *GetExclusions(Bead *bead);
  const exclusion_t *GetExclusions(Bead *bead) const;

  using iterator = std::list<exclusion_t *>::iterator;

  iterator begin() { return _exclusions.begin(); }
  iterator end() { return _exclusions.end(); }

  bool IsExcluded(Bead *bead1, Bead *bead2) const;

  template <typename iterable>
  void InsertExclusion(Bead *bead, iterable &excluded);

  void InsertExclusion(Bead *bead1, Bead *bead2);

  void RemoveExclusion(Bead *bead1, Bead *bead2);

 private:
  std::list<exclusion_t *> _exclusions;
  std::map<Bead *, exclusion_t *> _excl_by_bead;

  friend std::ostream &operator<<(std::ostream &out, ExclusionList &exl);
};

template <typename iterable>
inline void ExclusionList::Remove(iterable &l) {
  typename iterable::iterator i, j;

  for (i = l.begin(); i != l.end(); ++i) {
    for (j = i; j != l.end(); ++j) {
      RemoveExclusion(*i, *j);
    }
  }
}

template <typename iterable>
inline void ExclusionList::ExcludeList(iterable &l) {
  typename iterable::iterator i, j;

  for (i = l.begin(); i != l.end(); ++i) {
    for (j = i; j != l.end(); ++j) {
      InsertExclusion(*i, *j);
    }
  }
}

template <typename iterable>
inline void ExclusionList::InsertExclusion(Bead *beadA, iterable &l) {
  for (Bead *beadB : l) {
    Bead *bead1 = beadA;
    Bead *bead2 = beadB;
    if (bead2->getId() < bead1->getId()) {
      std::swap(bead1, bead2);
    }
    if (bead1 == bead2) {
      continue;
    }
    if (IsExcluded(bead1, bead2)) {
      continue;
    }

    exclusion_t *e;
    if ((e = GetExclusions(bead1)) == nullptr) {
      e = new exclusion_t;
      e->_atom = bead1;
      _exclusions.push_back(e);
      _excl_by_bead[bead1] = e;
    }
    e->_exclude.push_back(bead2);
  }
}

std::ostream &operator<<(std::ostream &out, ExclusionList &exl);

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_EXCLUSIONLIST_H
