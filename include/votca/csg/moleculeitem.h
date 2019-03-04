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

#ifndef _VOTCA_CSG_MOLECULEITEM_H
#define _VOTCA_CSG_MOLECULEITEM_H

#include <cassert>

namespace votca {
namespace csg {

class Molecule;

class MoleculeItem {
 public:
  virtual ~MoleculeItem() {}

  /**
   * Returns the molecule the pointer points at
   */
  Molecule *getMolecule() const {
    assert(_mol != nullptr);
    return _mol;
  }

  /**
   * stores a pointer to a molecule
   */
  void setMolecule(Molecule *mol) { _mol = mol; }

 protected:
  MoleculeItem(Molecule *mol) : _mol(mol) {}

  Molecule *_mol = nullptr;

  friend class BaseBead;
};
}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_MOLECULEITEM_H
