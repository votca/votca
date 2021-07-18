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

// Standard includes
#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

// Third party includes
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/property.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/csg/cgmoleculedef.h"
#include "votca/csg/interaction.h"
#include "votca/csg/map.h"
#include "votca/csg/topology.h"

namespace votca {
namespace csg {
class Molecule;
class Residue;
}  // namespace csg
}  // namespace votca

namespace votca {
namespace csg {

using namespace std;
using boost::lexical_cast;

CGMoleculeDef::~CGMoleculeDef() {
  for (beaddef_t *def : beads_) {
    delete def;
  }
  beads_.clear();
}

void CGMoleculeDef::Load(string filename) {
  options_.LoadFromXML(filename);
  // parse xml tree
  name_ = options_.get("cg_molecule.name").as<string>();
  ident_ = options_.get("cg_molecule.ident").as<string>();

  ParseTopology(options_.get("cg_molecule.topology"));
  ParseMapping(options_.get("cg_molecule.maps"));
}

void CGMoleculeDef::ParseTopology(tools::Property &options) {
  ParseBeads(options.get("cg_beads"));
  if (options.exists("cg_bonded")) {
    ParseBonded(options.get("cg_bonded"));
  }
}

void CGMoleculeDef::ParseBeads(tools::Property &options) {

  for (tools::Property *p : options.Select("cg_bead")) {
    beaddef_t *beaddef = new beaddef_t;
    beaddef->options_ = p;

    beaddef->name_ = p->get("name").as<string>();
    beaddef->type_ = p->get("type").as<string>();
    beaddef->mapping_ = p->get("mapping").as<string>();
    if (p->exists("symmetry")) {
      Index sym = p->get("symmetry").as<Index>();
      if (sym == 1) {
        beaddef->symmetry_ = Bead::spherical;
      } else if (sym == 3) {
        beaddef->symmetry_ = Bead::ellipsoidal;
      } else {
        throw std::runtime_error(
            "Only beads with spherical(1) or ellipsoidal(3) symmetry "
            "implemented.");
      }
    } else {
      beaddef->symmetry_ = Bead::spherical;
    }

    if (beads_by_name_.find(beaddef->name_) != beads_by_name_.end()) {
      throw std::runtime_error(string("bead name ") + beaddef->name_ +
                               " not unique in mapping");
    }
    beads_.push_back(beaddef);
    beads_by_name_[beaddef->name_] = beaddef;
  }
}

void CGMoleculeDef::ParseBonded(tools::Property &options) {
  bonded_ = options.Select("*");
}

void CGMoleculeDef::ParseMapping(tools::Property &options) {

  for (tools::Property *p : options.Select("map")) {
    maps_[p->get("name").as<string>()] = p;
  }
}
Molecule *CGMoleculeDef::CreateMolecule(Topology &top) {
  // add the residue names
  const Residue &res = top.CreateResidue(name_);
  Molecule *minfo = top.CreateMolecule(name_);

  // create the atoms
  for (auto &bead_def : beads_) {

    string type = bead_def->type_;
    if (!top.BeadTypeExist(type)) {
      top.RegisterBeadType(type);
    }
    Bead *bead = top.CreateBead(bead_def->symmetry_, bead_def->name_, type,
                                res.getId(), 0, 0);
    minfo->AddBead(bead, bead->getName());
  }

  // create the bonds
  map<string, string> had_iagroup;

  for (tools::Property *prop : bonded_) {
    std::list<Index> atoms;
    string iagroup = prop->get("name").as<string>();

    if (had_iagroup[iagroup] == "yes") {
      throw runtime_error(
          string("double occurence of interactions with name ") + iagroup);
    }
    had_iagroup[iagroup] = "yes";

    tools::Tokenizer tok(prop->get("beads").value(), " \n\t");
    for (auto &atom : tok) {
      Index i = minfo->getBeadIdByName(atom);
      if (i < 0) {
        throw runtime_error(
            string("error while trying to create bonded interaction, "
                   "bead " +
                   atom + " not found"));
      }

      atoms.push_back(i);
    }

    Index NrBeads = 1;
    if (prop->name() == "bond") {
      NrBeads = 2;
    } else if (prop->name() == "angle") {
      NrBeads = 3;
    } else if (prop->name() == "dihedral") {
      NrBeads = 4;
    }

    if ((atoms.size() % NrBeads) != 0) {
      throw runtime_error("Number of atoms in interaction '" +
                          prop->get("name").as<string>() +
                          "' is not a multiple of " +
                          lexical_cast<string>(NrBeads) + "! Missing beads?");
    }

    Index index = 0;
    while (!atoms.empty()) {
      Interaction *ic;

      if (prop->name() == "bond") {
        ic = new IBond(atoms);
      } else if (prop->name() == "angle") {
        ic = new IAngle(atoms);
      } else if (prop->name() == "dihedral") {
        ic = new IDihedral(atoms);
      } else {
        throw runtime_error("unknown bonded type in map: " + prop->name());
      }

      ic->setGroup(iagroup);
      ic->setIndex(index);
      ic->setMolecule(minfo->getId());
      top.AddBondedInteraction(ic);
      minfo->AddInteraction(ic);
      index++;
    }
  }
  return minfo;
}

Map CGMoleculeDef::CreateMap(const Molecule &in, Molecule &out) {
  if (out.BeadCount() != Index(beads_.size())) {
    throw runtime_error(
        "number of beads for cg molecule and mapping definition do "
        "not match, check your molecule naming.");
  }

  Map map(in, out);
  for (auto &bead : beads_) {

    Index iout = out.getBeadByName(bead->name_);
    if (iout < 0) {
      throw runtime_error(string("mapping error: reference molecule " +
                                 bead->name_ + " does not exist"));
    }

    tools::Property *mdef = getMapByName(bead->mapping_);
    if (!mdef) {
      throw runtime_error(string("mapping " + bead->mapping_ + " not found"));
    }

    /// TODO: change this to factory, do not hardcode!!
    BeadMap *bmap;
    switch (bead->symmetry_) {
      case 1:
        bmap = map.CreateBeadMap(BeadMapType::Spherical);
        break;
      case 3:
        bmap = map.CreateBeadMap(BeadMapType::Ellipsoidal);
        break;
      default:
        throw runtime_error(string("unknown symmetry in bead definition!"));
    }
    ////////////////////////////////////////////////////

    bmap->Initialize(&in, out.getBead(iout), (bead->options_), mdef);
  }
  return map;
}

CGMoleculeDef::beaddef_t *CGMoleculeDef::getBeadByName(const string &name) {
  map<string, beaddef_t *>::iterator iter = beads_by_name_.find(name);
  if (iter == beads_by_name_.end()) {
    std::cout << "cannot find: <" << name << "> in " << name_ << "\n";
    return nullptr;
  }
  return (*iter).second;
}

tools::Property *CGMoleculeDef::getMapByName(const string &name) {
  map<string, tools::Property *>::iterator iter = maps_.find(name);
  if (iter == maps_.end()) {
    std::cout << "cannot find map " << name << "\n";
    return nullptr;
  }
  return (*iter).second;
}

}  // namespace csg
}  // namespace votca
