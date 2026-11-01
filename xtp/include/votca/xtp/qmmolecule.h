/*
 *            Copyright 2009-2023 The VOTCA Development Team
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

#pragma once
#ifndef VOTCA_XTP_QMMOLECULE_H
#define VOTCA_XTP_QMMOLECULE_H

// Local VOTCA includes
#include "atomcontainer.h"
#include "qmatom.h"

namespace votca {
namespace xtp {

class QMMolecule : public AtomContainer<QMAtom> {
 public:
  QMMolecule(std::string name, Index id) : AtomContainer<QMAtom>(name, id) {};

  QMMolecule(CheckpointReader& r) : AtomContainer<QMAtom>(r) {};
  void LoadFromFile(std::string filename);

  void WriteXYZ(std::string filename, std::string header) const;

  void AddContainer(const AtomContainer<QMAtom>& container) {
    Index offset = atomlist_.size();
    type_ += "_" + container.getType();
    for (const auto& at : container) {
      // Copy the whole atom so that external-bond data and bonded partners
      // are kept (constructing a new QMAtom from index, element and
      // position dropped them). The atom ID and the bonded-partner IDs
      // are local to the added container and are shifted by offset; the
      // external-bond partner is a segment ID and needs no offset.
      QMAtom atom = at;
      atom.setID(at.getId() + offset);
      const Index* partners = at.getBondedPartnerIds();
      Index offset_partners[QMAtom::kMaxBondedPartners];
      for (Index i = 0; i < QMAtom::kMaxBondedPartners; i++) {
        offset_partners[i] = (partners[i] == -1) ? -1 : partners[i] + offset;
      }
      atom.setBondedPartnerIds(offset_partners);
      atomlist_.push_back(atom);
    }
    calcPos();
  }

  void ReorderAtomIDs() {
    Index id = 0;
    for (auto& at : atomlist_) {
      at.setID(id);
      id++;
    }
  }

  friend std::ostream& operator<<(std::ostream& out,
                                  const QMMolecule& container) {
    out << container.getId() << " " << container.getType() << "\n";
    for (const QMAtom& atom : container) {
      out << atom;
    }
    out << std::endl;
    return out;
  }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMMOLECULE_H
