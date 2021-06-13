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
#include <cassert>
#include <memory>
#include <regex>
#include <stdexcept>
#include <unordered_set>

// Third party includes
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/rangeparser.h>

// Local VOTCA includes
#include "votca/csg/boundarycondition.h"
#include "votca/csg/interaction.h"
#include "votca/csg/molecule.h"
#include "votca/csg/openbox.h"
#include "votca/csg/topology.h"

namespace votca {
namespace csg {

using namespace std;

bool is_digits(const std::string &str) {
  return str.find_first_not_of("0123456789") == std::string::npos;
}

Topology::~Topology() { Cleanup(); }

void Topology::Cleanup() {
  // cleanup beads
  beads_.clear();

  // cleanup molecules
  molecules_.clear();

  // cleanup residues
  residues_.clear();

  // cleanup interactions
  {
    for (InteractionContainer::iterator i = interactions_.begin();
         i < interactions_.end(); ++i) {
      delete (*i);
    }
    interactions_.clear();
  }
  // cleanup  bc_ object
  bc_ = std::make_unique<OpenBox>();
}

/// \todo implement checking, only used in xml topology reader
void Topology::CreateMoleculesByRange(string name, Index first, Index nbeads,
                                      Index nmolecules) {
  Molecule *mol = CreateMolecule(name);
  Index beadcount = 0;
  Index res_offset = 0;

  for (auto &bead_ : beads_) {
    // xml numbering starts with 1
    if (--first > 0) {
      continue;
    }
    // This is not 100% correct, but let's assume for now that the resnr do
    // increase
    if (beadcount == 0) {
      res_offset = bead_.getResnr();
    }
    stringstream bname;
    bname << bead_.getResnr() - res_offset + 1 << ":"
          << getResidue(bead_.getResnr()).getName() << ":" << bead_.getName();
    mol->AddBead(&bead_, bname.str());
    if (++beadcount == nbeads) {
      if (--nmolecules <= 0) {
        break;
      }
      mol = CreateMolecule(name);
      beadcount = 0;
    }
  }
}

/// \todo clean up CreateMoleculesByResidue!
void Topology::CreateMoleculesByResidue() {
  // first create a molecule for each residue
  for (const auto &residue_ : residues_) {
    CreateMolecule(residue_.getName());
  }

  // add the beads to the corresponding molecules based on their resid
  for (auto &bead_ : beads_) {

    MoleculeByIndex(bead_.getResnr())
        ->AddBead(&bead_, string("1:TRI:") + bead_.getName());
  }

  /// \todo sort beads in molecules that all beads are stored in the same order.
  /// This is needed for the mapping!
}

void Topology::CreateOneBigMolecule(string name) {
  Molecule *mi = CreateMolecule(name);

  for (auto &bead_ : beads_) {
    stringstream n("");
    n << bead_.getResnr() + 1 << ":" << residues_[bead_.getResnr()].getName()
      << ":" << bead_.getName();
    mi->AddBead(&bead_, n.str());
  }
}

void Topology::Add(Topology *top) {

  Index res0 = ResidueCount();

  for (const auto &bi : top->beads_) {
    string type = bi.getType();
    CreateBead(bi.getSymmetry(), bi.getName(), type, bi.getResnr() + res0,
               bi.getMass(), bi.getQ());
  }

  for (const auto &residue_ : top->residues_) {
    CreateResidue(residue_.getName());
  }

  // \todo beadnames in molecules!!
  for (auto &molecule_ : top->molecules_) {
    Molecule *mi = CreateMolecule(molecule_.getName());
    for (Index i = 0; i < mi->BeadCount(); i++) {
      mi->AddBead(mi->getBead(i), "invalid");
    }
  }
}

void Topology::CopyTopologyData(Topology *top) {

  bc_->setBox(top->getBox());
  time_ = top->time_;
  step_ = top->step_;

  // cleanup old data
  Cleanup();

  // copy all residues
  for (const auto &residue_ : top->residues_) {
    CreateResidue(residue_.getName());
  }

  // create all beads
  for (const auto &bi : top->beads_) {
    string type = bi.getType();
    CreateBead(bi.getSymmetry(), bi.getName(), type, bi.getResnr(),
               bi.getMass(), bi.getQ());
  }

  // copy all molecules
  for (const auto &molecule_ : top->molecules_) {
    Molecule *mi = CreateMolecule(molecule_.getName());
    for (Index i = 0; i < molecule_.BeadCount(); i++) {
      Index beadid = molecule_.getBead(i)->getId();
      mi->AddBead(&beads_[beadid], molecule_.getBeadName(i));
    }
  }
}

Index Topology::getBeadTypeId(string type) const {
  assert(beadtypes_.count(type));
  return beadtypes_.at(type);
}

void Topology::RenameMolecules(string range, string name) {
  tools::RangeParser rp;
  rp.Parse(range);
  for (Index i : rp) {
    if (i > Index(molecules_.size())) {
      throw runtime_error(
          string("RenameMolecules: num molecules smaller than"));
    }
    getMolecule(i - 1)->setName(name);
  }
}

void Topology::RenameBeadType(string name, string newname) {

  for (auto &bead : beads_) {
    string type = bead.getType();
    if (tools::wildcmp(name, type)) {
      bead.setType(newname);
    }
  }
}

void Topology::SetBeadTypeMass(string name, double value) {
  for (auto &bead : beads_) {
    string type = bead.getType();
    if (tools::wildcmp(name, type)) {
      bead.setMass(value);
    }
  }
}

void Topology::CheckMoleculeNaming(void) {
  map<string, Index> nbeads;

  for (const auto &mol : molecules_) {
    map<string, Index>::iterator entry = nbeads.find(mol.getName());
    if (entry != nbeads.end()) {
      if (entry->second != mol.BeadCount()) {
        throw runtime_error(
            "There are molecules which have the same name but different number "
            "of bead "
            "please check the section manual topology handling in the votca "
            "manual");
      }
      continue;
    }
    nbeads[mol.getName()] = mol.BeadCount();
  }
}

void Topology::AddBondedInteraction(Interaction *ic) {
  map<string, Index>::iterator iter;
  iter = interaction_groups_.find(ic->getGroup());
  if (iter != interaction_groups_.end()) {
    ic->setGroupId((*iter).second);
  } else {
    Index i = interaction_groups_.size();
    interaction_groups_[ic->getGroup()] = i;
    ic->setGroupId(i);
  }
  interactions_.push_back(ic);
  interactions_by_group_[ic->getGroup()].push_back(ic);
}

std::vector<Interaction *> Topology::InteractionsInGroup(const string &group) {
  map<string, vector<Interaction *>>::iterator iter =
      interactions_by_group_.find(group);
  if (iter == interactions_by_group_.end()) {
    return vector<Interaction *>();
  }
  return iter->second;
}

bool Topology::BeadTypeExist(string type) const {
  return beadtypes_.count(type);
}

void Topology::RegisterBeadType(string type) {
  unordered_set<Index> ids;
  for (pair<const string, Index> type_and_id : beadtypes_) {
    ids.insert(type_and_id.second);
  }

  Index id = 0;
  // If the type is also a number use it as the id as well provided it is not
  // already taken
  if (is_digits(type)) {
    id = boost::lexical_cast<Index>(type);
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
  return bc_->BCShortestConnection(r_i, r_j);
}

Eigen::Vector3d Topology::getDist(Index bead1, Index bead2) const {
  return BCShortestConnection(getBead(bead1)->getPos(),
                              getBead(bead2)->getPos());
}

double Topology::BoxVolume() const { return bc_->BoxVolume(); }

void Topology::RebuildExclusions() { exclusions_.CreateExclusions(this); }

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
  return bc_->getShortestBoxDimension();
}

}  // namespace csg
}  // namespace votca
