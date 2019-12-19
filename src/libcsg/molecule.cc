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

#include <iostream>
#include <votca/csg/molecule.h>

namespace votca {
namespace csg {

using namespace std;

void Molecule::AddBead(Bead *bead, const string &name) {
  _beads.push_back(bead);
  _bead_names.push_back(name);
  _beadmap[name] = _beads.size() - 1;

  bead->setMoleculeId(_id);
}

Index Molecule::getBeadByName(const string &name) {
  map<string, Index>::iterator iter = _beadmap.find(name);
  if (iter == _beadmap.end()) {
    std::cout << "cannot find: <" << name << "> in " << _name << "\n";
    return -1;
  }
  return _beadmap[name];
}

}  // namespace csg
}  // namespace votca
