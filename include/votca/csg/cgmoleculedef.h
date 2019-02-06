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

#ifndef _VOTCA_CSG_CGMOLECULEDEF_H
#define _VOTCA_CSG_CGMOLECULEDEF_H

#include <list>
#include <map>
#include <string>
#include <vector>

#include "exclusionlist.h"
#include "map.h"
#include "molecule.h"
#include <votca/tools/property.h>
#include <votca/tools/types.h>

namespace votca {
namespace csg {
using namespace votca::tools;

/**
    \brief definition of a coarse grained molecule

    This class is to define a coarse grained molecule, which includes the
   topology, mapping, ...

    \todo clean up this class, do the bonded interactions right!!!!
    \todo check for consistency of xml file, seperate xml parser and class!!
*/
class CGMoleculeDef {
 public:
  CGMoleculeDef() {}
  ~CGMoleculeDef();

  Molecule *CreateMolecule(Topology &top);
  Map *CreateMap(Molecule &in, Molecule &out);

  void Load(std::string filename);

  const std::string &getName() { return _name; }
  const std::string &getIdent() { return _ident; }

 private:
  Property _options;

  struct beaddef_t {
    std::string _name;
    std::string _type;
    byte_t _symmetry;
    std::string _mapping;
    std::vector<std::string> _subbeads;
    Property *_options;
  };

  // name of the coarse grained molecule
  std::string _name;
  // name of the molecule to coarse grain
  std::string _ident;

  // beads of the cg molecule
  std::vector<beaddef_t *> _beads;
  std::map<std::string, beaddef_t *> _beads_by_name;

  // mapping schemes
  std::map<std::string, Property *> _maps;

  std::list<Property *> _bonded;

  void ParseTopology(Property &options);
  void ParseBeads(Property &options);
  void ParseBonded(Property &options);
  void ParseMapping(Property &options);

  beaddef_t *getBeadByName(const std::string &name);
  Property *getMapByName(const std::string &name);
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_CGMOLECULEDEF_H */
