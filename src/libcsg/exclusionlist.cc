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

  for (auto &exclusion_ : exclusions_) {
    delete exclusion_;
  }
  exclusions_.clear();
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
      excl_by_bead_.find(bead);
  if (iter == excl_by_bead_.end()) {
    return nullptr;
  }

  return (*iter).second;
}

ExclusionList::exclusion_t *ExclusionList::GetExclusions(Bead *bead) {
  std::map<Bead *, exclusion_t *>::iterator iter = excl_by_bead_.find(bead);
  if (iter == excl_by_bead_.end()) {
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
    if (find(excl->exclude_.begin(), excl->exclude_.end(), bead2) !=
        excl->exclude_.end()) {
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
    e->atom_ = bead1;
    exclusions_.push_back(e);
    excl_by_bead_[bead1] = e;
  }
  e->exclude_.push_back(bead2);
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
      std::find_if(exclusions_.begin(), exclusions_.end(),
                   [&bead1](exclusion_t *e) { return e->atom_ == bead1; });

  if (ex == exclusions_.end()) {
    return;
  }

  (*ex)->exclude_.remove(bead2);
  if ((*ex)->exclude_.empty()) {
    (*ex) = nullptr;
    exclusions_.erase(ex);
  }
  exclusions_.remove(nullptr);
}

bool compareAtomIdiExclusionList(const ExclusionList::exclusion_t *a,
                                 const ExclusionList::exclusion_t *b) {
  return a->atom_->getId() < b->atom_->getId();
}

bool compareAtomIdBeadList(const Bead *a, const Bead *b) {
  return a->getId() < b->getId();
}

std::ostream &operator<<(std::ostream &out, ExclusionList &exl) {
  exl.exclusions_.sort(compareAtomIdiExclusionList);

  for (auto &exclusion_ : exl.exclusions_) {
    exclusion_->exclude_.sort(compareAtomIdBeadList);
    out << (Index)(exclusion_->atom_->getId()) + 1;
    for (Bead *bead : exclusion_->exclude_) {
      out << " " << (bead->getId() + 1);
    }
    out << endl;
  }
  return out;
}

}  // namespace csg
}  // namespace votca
