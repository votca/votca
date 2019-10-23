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

#ifndef _VOTCA_CSG_RESIDUE_H
#define _VOTCA_CSG_RESIDUE_H

#include "topologyitem.h"
#include <string>

namespace votca {
namespace csg {

/**
    \brief class for a residue

    The Residue class describes a residue. When reading in the atoms, all beads
   belong to a residue. Later on, the molecules can be organized into molecules
   based on their residue.

*/
class Residue : public TopologyItem {
 public:
  /// get the name of the residue
  const std::string &getName();

  /// get the name of the residue
  const long int &getId() const { return _id; }

 private:
  long int _id;
  std::string _name;

 private:
  /// constructor
  Residue(Topology *parent, long int id, const std::string &name)
      : TopologyItem(parent), _id(id), _name(name) {}
  friend class Topology;
};

inline const std::string &Residue::getName() { return _name; }

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_RESIDUE_H */
