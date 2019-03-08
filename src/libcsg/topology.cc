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

#include <boost/lexical_cast.hpp>
#include <cassert>
#include <regex>
#include <stdexcept>
#include <unordered_set>
#include <votca/csg/boundarycondition.h>
#include <votca/csg/interaction.h>
#include <votca/csg/molecule.h>
#include <votca/csg/openbox.h>
#include <votca/csg/topology.h>
#include <votca/tools/rangeparser.h>

namespace votca {
namespace csg {

using namespace std;

bool is_digits(const std::string &str) {
  return str.find_first_not_of("0123456789") == std::string::npos;
}

Topology::~Topology() {
  Cleanup();
  if (_bc) delete (_bc);
  _bc = nullptr;
}

void Topology::Cleanup() {
  // cleanup beads
  {
    BeadContainer::iterator i;
    for (i = _beads.begin(); i < _beads.end(); ++i) delete *i;
    _beads.clear();
  }
  // cleanup molecules
  {
    MoleculeContainer::iterator i;
    for (i = _molecules.begin(); i < _molecules.end(); ++i) delete *i;
    _molecules.clear();
  }
  // cleanup residues
  {
    ResidueContainer::iterator i;
    for (i = _residues.begin(); i < _residues.end(); ++i) delete (*i);
    _residues.clear();
  }
  // cleanup interactions
  {
    InteractionContainer::iterator i;
    for (i = _interactions.begin(); i < _interactions.end(); ++i) delete (*i);
    _interactions.clear();
  }
  // cleanup _bc object
  if (_bc) delete (_bc);
  _bc = new OpenBox();
}

/// \todo implement checking, only used in xml topology reader
void Topology::CreateMoleculesByRange(string name, int first, int nbeads,
                                      int nmolecules) {
  Molecule *mol        = CreateMolecule(name);
  int       beadcount  = 0;
  int       res_offset = 0;

  BeadContainer::iterator bead;
  for (bead = _beads.begin(); bead != _beads.end(); ++bead) {
    // xml numbering starts with 1
    if (--first > 0) continue;
    // This is not 100% correct, but let's assume for now that the resnr do
    // increase
    if (beadcount == 0) {
      res_offset = (*bead)->getResnr();
    }
    stringstream bname;
    bname << (*bead)->getResnr() - res_offset + 1 << ":"
          << getResidue((*bead)->getResnr())->getName() << ":"
          << (*bead)->getName();
    mol->AddBead((*bead), bname.str());
    if (++beadcount == nbeads) {
      if (--nmolecules <= 0) break;
      mol       = CreateMolecule(name);
      beadcount = 0;
    }
  }
}

/// \todo clean up CreateMoleculesByResidue!
void Topology::CreateMoleculesByResidue() {
  // first create a molecule for each residue
  ResidueContainer::iterator res;
  for (res = _residues.begin(); res != _residues.end(); ++res) {
    CreateMolecule((*res)->getName());
  }

  // add the beads to the corresponding molecules based on their resid
  BeadContainer::iterator bead;
  for (bead = _beads.begin(); bead != _beads.end(); ++bead) {

    MoleculeByIndex((*bead)->getResnr())
        ->AddBead((*bead), string("1:TRI:") + (*bead)->getName());
  }

  /// \todo sort beads in molecules that all beads are stored in the same order.
  /// This is needed for the mapping!
}

void Topology::CreateOneBigMolecule(string name) {
  Molecule *mi = CreateMolecule(name);

  BeadContainer::iterator bead;

  for (bead = _beads.begin(); bead != _beads.end(); ++bead) {
    stringstream n("");
    n << (*bead)->getResnr() + 1 << ":"
      << _residues[(*bead)->getResnr()]->getName() << ":" << (*bead)->getName();
    mi->AddBead((*bead), n.str());
  }
}

void Topology::Add(Topology *top) {
  BeadContainer::iterator     bead;
  ResidueContainer::iterator  res;
  MoleculeContainer::iterator mol;

  int res0 = ResidueCount();

  for (bead = top->_beads.begin(); bead != top->_beads.end(); ++bead) {
    Bead * bi   = *bead;
    string type = bi->getType();
    CreateBead(bi->getSymmetry(), bi->getName(), type, bi->getResnr() + res0,
               bi->getMass(), bi->getQ());
  }

  for (res = top->_residues.begin(); res != top->_residues.end(); ++res) {
    CreateResidue((*res)->getName());
  }

  // \todo beadnames in molecules!!
  for (mol = top->_molecules.begin(); mol != top->_molecules.end(); ++mol) {
    Molecule *mi = CreateMolecule((*mol)->getName());
    for (int i = 0; i < mi->BeadCount(); i++) {
      mi->AddBead(mi->getBead(i), "invalid");
    }
  }
}

void Topology::CopyTopologyData(Topology *top) {
  BeadContainer::iterator     it_bead;
  ResidueContainer::iterator  it_res;
  MoleculeContainer::iterator it_mol;

  _bc->setBox(top->getBox());
  _time = top->_time;
  _step = top->_step;

  // cleanup old data
  Cleanup();

  // copy all residues
  for (it_res = top->_residues.begin(); it_res != top->_residues.end();
       ++it_res) {
    CreateResidue((*it_res)->getName());
  }

  // create all beads
  for (it_bead = top->_beads.begin(); it_bead != top->_beads.end(); ++it_bead) {
    Bead * bi   = *it_bead;
    string type = bi->getType();
    CreateBead(bi->getSymmetry(), bi->getName(), type, bi->getResnr(),
               bi->getMass(), bi->getQ());
  }

  // copy all molecules
  for (it_mol = top->_molecules.begin(); it_mol != top->_molecules.end();
       ++it_mol) {
    Molecule *mi = CreateMolecule((*it_mol)->getName());
    for (int i = 0; i < (*it_mol)->BeadCount(); i++) {
      int beadid = (*it_mol)->getBead(i)->getId();
      mi->AddBead(_beads[beadid], (*it_mol)->getBeadName(i));
    }
  }
}

int Topology::getBeadTypeId(string type) const {
  assert(beadtypes_.count(type));
  return beadtypes_.at(type);
}

void Topology::RenameMolecules(string range, string name) {
  RangeParser           rp;
  RangeParser::iterator i;

  rp.Parse(range);
  for (i = rp.begin(); i != rp.end(); ++i) {
    if ((unsigned int)*i > _molecules.size())
      throw runtime_error(
          string("RenameMolecules: num molecules smaller than"));
    getMolecule(*i - 1)->setName(name);
  }
}

void Topology::RenameBeadType(string name, string newname) {
  BeadContainer::iterator bead;
  for (bead = _beads.begin(); bead != _beads.end(); ++bead) {
    string type = (*bead)->getType();
    if (wildcmp(name.c_str(), type.c_str())) {
      (*bead)->setType(newname);
    }
  }
}

void Topology::SetBeadTypeMass(string name, double value) {
  BeadContainer::iterator bead;
  for (bead = _beads.begin(); bead != _beads.end(); ++bead) {
    string type = (*bead)->getType();
    if (wildcmp(name.c_str(), type.c_str())) {
      (*bead)->setMass(value);
    }
  }
}

void Topology::CheckMoleculeNaming(void) {
  map<string, int> nbeads;

  for (MoleculeContainer::iterator iter = _molecules.begin();
       iter != _molecules.end(); ++iter) {
    map<string, int>::iterator entry = nbeads.find((*iter)->getName());
    if (entry != nbeads.end()) {
      if (entry->second != (*iter)->BeadCount())
        throw runtime_error(
            "There are molecules which have the same name but different number "
            "of bead "
            "please check the section manual topology handling in the votca "
            "manual");
      continue;
    }
    nbeads[(*iter)->getName()] = (*iter)->BeadCount();
  }
}

void Topology::AddBondedInteraction(Interaction *ic) {
  map<string, int>::iterator iter;
  iter = _interaction_groups.find(ic->getGroup());
  if (iter != _interaction_groups.end())
    ic->setGroupId((*iter).second);
  else {
    int i                               = _interaction_groups.size();
    _interaction_groups[ic->getGroup()] = i;
    ic->setGroupId(i);
  }
  _interactions.push_back(ic);
  _interactions_by_group[ic->getGroup()].push_back(ic);
}

std::list<Interaction *> Topology::InteractionsInGroup(const string &group) {
  map<string, list<Interaction *>>::iterator iter;
  iter = _interactions_by_group.find(group);
  if (iter == _interactions_by_group.end()) return list<Interaction *>();
  return iter->second;
}

bool Topology::BeadTypeExist(string type) const {
  return beadtypes_.count(type);
}

void Topology::RegisterBeadType(string type) {
  unordered_set<int> ids;
  for (pair<const string, int> type_and_id : beadtypes_) {
    ids.insert(type_and_id.second);
  }

  int id = 0;
  // If the type is also a number use it as the id as well provided it is not
  // already taken
  if (is_digits(type)) {
    id = boost::lexical_cast<int>(type);
    assert(!ids.count(id) &&
           "The type passed in is a number and has already"
           " been registered. It is likely that you are passing in numbers as "
           "bead types as well as strings, choose one or the other do not mix "
           "between using numbers and strings ");
  }

  while (ids.count(id)) {
    ++id;
  }
  beadtypes_[type] = id;
}

Eigen::Vector3d Topology::BCShortestConnection(
    const Eigen::Vector3d &r_i, const Eigen::Vector3d &r_j) const {
  return _bc->BCShortestConnection(r_i, r_j);
}

Eigen::Vector3d Topology::getDist(int bead1, int bead2) const {
  return BCShortestConnection(getBead(bead1)->getPos(),
                              getBead(bead2)->getPos());
}

double Topology::BoxVolume() const { return _bc->BoxVolume(); }

void Topology::RebuildExclusions() { _exclusions.CreateExclusions(this); }

BoundaryCondition::eBoxtype Topology::autoDetectBoxType(
    const Eigen::Matrix3d &box) const {
  // set the box type to OpenBox in case "box" is the zero matrix,
  // to OrthorhombicBox in case "box" is a diagonal matrix,
  // or to TriclinicBox otherwise
  if (box.isApproxToConstant(0)) {
    return BoundaryCondition::typeOpen;
  } else if ((box - Eigen::Matrix3d(box.diagonal().asDiagonal()))
                 .isApproxToConstant(0)) {
    return BoundaryCondition::typeOrthorhombic;
  } else {
    return BoundaryCondition::typeTriclinic;
  }
  return BoundaryCondition::typeOpen;
}

double Topology::ShortestBoxSize() const {
  Eigen::Vector3d box_a = getBox().col(0);
  Eigen::Vector3d box_b = getBox().col(1);
  Eigen::Vector3d box_c = getBox().col(2);

  // create plane normals
  Eigen::Vector3d norm_a = box_b.cross(box_c);
  Eigen::Vector3d norm_b = box_c.cross(box_a);
  Eigen::Vector3d norm_c = box_a.cross(box_b);

  norm_a.normalize();
  norm_b.normalize();
  norm_c.normalize();

  double la = box_a.dot(norm_a);
  double lb = box_b.dot(norm_b);
  double lc = box_c.dot(norm_c);

  return std::min(la, std::min(lb, lc));
}

}  // namespace csg
}  // namespace votca
