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

// Standard includes
#include <algorithm>

// Local VOTCA includes
#include "votca/csg/exclusionlist.h"
#include "votca/csg/topology.h"

namespace votca {
namespace csg {

using namespace std;

void ExclusionList::Clear(void) {

  for (auto &_exclusion : _exclusions) {
    delete _exclusion;
  }
  _exclusions.clear();
}

void ExclusionList::CreateExclusions(Topology *top) {
  InteractionContainer &ic = top->BondedInteractions();

  for (auto &ia : ic) {
    Index beads_in_int = ia->BeadCount();
    std::vector<Bead *> l;

    for (Index ibead = 0; ibead < beads_in_int; ibead++) {
      Index ii = ia->getBeadId(ibead);
      l.push_back(top->getBead(ii));
    }
    ExcludeList(l);
  }
}

const ExclusionList::exclusion_t *ExclusionList::GetExclusions(
    Bead *bead) const {
  std::map<Bead *, exclusion_t *>::const_iterator iter =
      _excl_by_bead.find(bead);
  if (iter == _excl_by_bead.end()) {
    return nullptr;
  }

  return (*iter).second;
}

ExclusionList::exclusion_t *ExclusionList::GetExclusions(Bead *bead) {
  std::map<Bead *, exclusion_t *>::iterator iter = _excl_by_bead.find(bead);
  if (iter == _excl_by_bead.end()) {
    return nullptr;
  }

  return (*iter).second;
}

bool ExclusionList::IsExcluded(Bead *bead1, Bead *bead2) const {
  if (bead1->getMoleculeId() != bead2->getMoleculeId()) {
    return false;
  }
  if (bead2->getId() < bead1->getId()) {
    swap(bead1, bead2);
  }

  const exclusion_t *excl = GetExclusions(bead1);
  if (excl != nullptr) {
    if (find(excl->_exclude.begin(), excl->_exclude.end(), bead2) !=
        excl->_exclude.end()) {
      return true;
    }
  }
  return false;
}

void ExclusionList::InsertExclusion(Bead *bead1, Bead *bead2) {
  if (bead2->getId() < bead1->getId()) {
    std::swap(bead1, bead2);
  }

  if (bead1 == bead2) {
    return;
  }

  if (IsExcluded(bead1, bead2)) {
    return;
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

void ExclusionList::RemoveExclusion(Bead *bead1, Bead *bead2) {
  if (bead2->getId() < bead1->getId()) {
    std::swap(bead1, bead2);
  }

  if (bead1 == bead2) {
    return;
  }

  if (!IsExcluded(bead1, bead2)) {
    return;
  }

  std::list<exclusion_t *>::iterator ex =
      std::find_if(_exclusions.begin(), _exclusions.end(),
                   [&bead1](exclusion_t *e) { return e->_atom == bead1; });

  if (ex == _exclusions.end()) {
    return;
  }

  (*ex)->_exclude.remove(bead2);
  if ((*ex)->_exclude.empty()) {
    (*ex) = nullptr;
    _exclusions.erase(ex);
  }
  _exclusions.remove(nullptr);
}

bool compareAtomIdiExclusionList(const ExclusionList::exclusion_t *a,
                                 const ExclusionList::exclusion_t *b) {
  return a->_atom->getId() < b->_atom->getId();
}

bool compareAtomIdBeadList(const Bead *a, const Bead *b) {
  return a->getId() < b->getId();
}

std::ostream &operator<<(std::ostream &out, ExclusionList &exl) {
  exl._exclusions.sort(compareAtomIdiExclusionList);

  for (auto &_exclusion : exl._exclusions) {
    _exclusion->_exclude.sort(compareAtomIdBeadList);
    out << (Index)(_exclusion->_atom->getId()) + 1;
    for (Bead *bead : _exclusion->_exclude) {
      out << " " << (bead->getId() + 1);
    }
    out << endl;
  }
  return out;
}

}  // namespace csg
}  // namespace votca
