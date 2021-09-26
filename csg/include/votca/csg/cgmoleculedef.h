/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_CSG_CGMOLECULEDEF_H
#define VOTCA_CSG_CGMOLECULEDEF_H

// Standard includes
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

// VOTCA includes
#include <votca/tools/property.h>
#include <votca/tools/types.h>

// Local VOTCA includes
#include "exclusionlist.h"
#include "map.h"
#include "molecule.h"

namespace votca {
namespace csg {

/**
    \brief definition of a coarse grained molecule

    This class is to define a coarse grained molecule, which includes the
   topology, mapping, ...

    \todo clean up this class, do the bonded interactions right!!!!
    \todo check for consistency of xml file, seperate xml parser and class!!
*/
class CGMoleculeDef {
 public:
  CGMoleculeDef() = default;
  ~CGMoleculeDef();

  Molecule *CreateMolecule(Topology &top);
  Map CreateMap(const Molecule &in, Molecule &out);

  void Load(std::string filename);

  const std::string &getName() { return name_; }
  const std::string &getIdent() { return ident_; }

 private:
  tools::Property options_;

  struct beaddef_t {
    std::string name_;
    std::string type_;
    Bead::Symmetry symmetry_;
    std::string mapping_;
    std::vector<std::string> subbeads_;
    tools::Property *options_;
  };

  // name of the coarse grained molecule
  std::string name_;
  // name of the molecule to coarse grain
  std::string ident_;

  // beads of the cg molecule
  std::vector<beaddef_t *> beads_;
  std::map<std::string, beaddef_t *> beads_by_name_;

  // mapping schemes
  std::map<std::string, tools::Property *> maps_;

  std::vector<tools::Property *> bonded_;

  void ParseTopology(tools::Property &options);
  void ParseBeads(tools::Property &options);
  void ParseBonded(tools::Property &options);
  void ParseMapping(tools::Property &options);

  beaddef_t *getBeadByName(const std::string &name);
  tools::Property *getMapByName(const std::string &name);
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CGMOLECULEDEF_H
