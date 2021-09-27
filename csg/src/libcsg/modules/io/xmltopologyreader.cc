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
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

// Third party includes
#include <boost/lexical_cast.hpp>

// Local private VOTCA includes
#include "xmltopologyreader.h"

namespace votca {
namespace csg {

using namespace std;

bool XMLTopologyReader::ReadTopology(string filename, Topology &top) {
  top_ = &top;

  tools::Property options;
  options.LoadFromXML(filename);
  ParseRoot(options.get("topology"));

  top_->RebuildExclusions();
  return true;
}

void XMLTopologyReader::ReadTopolFile(string file) {
  std::unique_ptr<TopologyReader> reader = TopReaderFactory().Create(file);
  if (!reader) {
    throw runtime_error(file + ": unknown topology format");
  }

  reader->ReadTopology(file, *top_);
  // Clean XML molecules and beads.
}

void XMLTopologyReader::ParseRoot(tools::Property &property) {
  has_base_topology_ = false;
  if (property.hasAttribute("base")) {
    ReadTopolFile(property.getAttribute<string>("base"));
    has_base_topology_ = true;
  }

  // Iterate over keys at first level.
  for (auto &it : property) {
    if (it.name() == "h5md_particle_group") {
      top_->setParticleGroup(it.getAttribute<string>("name"));
    } else if (it.name() == "molecules") {
      mol_index_ = 1;
      bead_index_ = 0;
      ParseMolecules(it);
    } else if (it.name() == "bonded") {
      ParseBonded(it);
    } else if (it.name() == "box") {
      ParseBox(it);
    } else if (it.name() == "beadtypes") {
      ParseBeadTypes(it);
    } else {
      throw runtime_error("unknown tag: topology." + it.name());
    }
  }
}

void XMLTopologyReader::ParseBox(tools::Property &p) {
  Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
  m(0, 0) = p.getAttribute<double>("xx");
  m(1, 1) = p.getAttribute<double>("yy");
  m(2, 2) = p.getAttribute<double>("zz");
  top_->setBox(m);
}

void XMLTopologyReader::ParseMolecules(tools::Property &p) {
  for (auto &it : p) {
    if (it.name() == "clear") {
      top_->ClearMoleculeList();
    } else if (it.name() == "rename") {
      string molname = it.getAttribute<string>("name");
      string range = it.getAttribute<string>("range");
      top_->RenameMolecules(range, molname);
    } else if (it.name() == "define" || it.name() == "molecule") {
      string molname = it.getAttribute<string>("name");
      Index first = 0;
      if (it.name() == "define") {
        first = it.getAttribute<Index>("first");
      }
      Index nbeads = it.getAttribute<Index>("nbeads");
      Index nmols = it.getAttribute<Index>("nmols");
      if (it.name() == "define" && first < 1) {
        throw std::runtime_error(
            "Attribute first is suppose to be > 0, but found " +
            boost::lexical_cast<string>(it.getAttribute<string>("first")));
      }
      if (nbeads < 1) {
        throw std::runtime_error(
            "Attribute nbeads is suppose to be > 0, but found " +
            boost::lexical_cast<string>(it.getAttribute<string>("nbeads")));
      }
      if (nmols < 1) {
        throw std::runtime_error(
            "Attribute nmols is suppose to be > 0, but found " +
            boost::lexical_cast<string>(it.getAttribute<string>("nmols")));
      }
      if (it.name() == "define") {
        top_->CreateMoleculesByRange(molname, first, nbeads, nmols);
      } else {
        if (has_base_topology_) {
          throw std::runtime_error(
              "The defined list of beads only works for pure xml topology, "
              "without 'base' attribute.");
        }
        ParseMolecule(it, molname, nmols);
      }
    }
  }
}

void XMLTopologyReader::ParseMolecule(tools::Property &p, string molname,
                                      Index nmols) {
  vector<XMLBead *> xmlBeads;
  vector<Index> xmlResidues;
  for (auto &it : p) {
    if (it.name() == "bead") {
      string atname = it.getAttribute<string>("name");
      string attype = it.getAttribute<string>("type");
      double atmass, atq;
      Index resid;
      try {
        atmass = it.getAttribute<double>("mass");
      } catch (runtime_error &) {
        atmass = 1.0;
      }
      try {
        atq = it.getAttribute<double>("q");
      } catch (runtime_error &) {
        atq = 0.0;
      }
      try {
        resid = it.getAttribute<Index>("resid");
        if (resid <= 0) {
          throw std::invalid_argument(
              "Residue count for beads in topology.molecules.molecule has to "
              "be greater than zero");
        }
      } catch (runtime_error &) {
        resid = -1;
      }
      if (!xmlResidues.empty()) {
        if (xmlResidues.back() != resid && xmlResidues.back() != resid - 1) {
          throw std::invalid_argument(
              "Residue count for beads in topology.molecules.molecule does not "
              "increase monotonically in steps of 1");
        }
        if (xmlResidues.back() != -1 && resid == -1) {
          throw std::invalid_argument(
              "Residue count for beads in topology.molecules.molecule has to "
              "be declared for all beads or for none");
        }
      }
      xmlResidues.push_back(resid);
      XMLBead *xmlBead = new XMLBead(atname, attype, atmass, atq);
      xmlBeads.push_back(xmlBead);
    } else {
      throw std::runtime_error(
          "Wrong element under topology.molecules.molecule: " + it.name());
    }
  }
  if (xmlResidues.size() != xmlBeads.size()) {
    throw std::runtime_error(
        "Number of elements in bead-vector and residue-vector are not "
        "identical");
  }
  // Create molecule in topology. Replicate data.
  Index resnr = top_->ResidueCount();
  if (!xmlResidues.empty()) {
    if (xmlResidues.front() != resnr + 1 && xmlResidues.front() != -1) {
      throw std::runtime_error(
          "Residue count for beads in topology.molecules.molecule has to be "
          "greater than the number of residues already in the topology");
    }
  }
  for (Index mn = 0; mn < nmols; mn++) {
    Molecule *mi = top_->CreateMolecule(molname);
    XMLMolecule *xmlMolecule = new XMLMolecule(molname, nmols);
    xmlMolecule->pid = mi->getId();
    xmlMolecule->mi = mi;
    molecules_.insert(make_pair(molname, xmlMolecule));
    vector<Index>::iterator resit = xmlResidues.begin();
    for (vector<XMLBead *>::iterator itb = xmlBeads.begin();
         itb != xmlBeads.end(); ++itb, ++resit) {
      stringstream bname;
      XMLBead &b = **itb;
      if (*resit != -1) {
        if (top_->ResidueCount() < *resit) {
          resnr = *resit - 1;
          top_->CreateResidue(molname, resnr);
        }
      } else {
        top_->CreateResidue(molname, resnr);
      }

      if (!top_->BeadTypeExist(b.type)) {
        top_->RegisterBeadType(b.type);
      }
      Bead *bead =
          top_->CreateBead(Bead::spherical, b.name, b.type, resnr, b.mass, b.q);
      bname << mol_index_ << ":" << molname << ":" << b.name;
      mi->AddBead(bead, bname.str());

      // Data for bonded terms.
      XMLBead *b_rep = new XMLBead(b);
      b_rep->pid = bead_index_;
      if (xmlMolecule->name2beads.count(b.name) != 0) {
        throw std::runtime_error("Atom " + b.name + " in molecule " + molname +
                                 " already exists.");
      }
      xmlMolecule->name2beads.insert(make_pair(b.name, b_rep));
      bead_index_++;
    }
    resnr++;
  }

  // clean up
  for (auto &xmlBead : xmlBeads) {
    delete xmlBead;
  }
}

void XMLTopologyReader::ParseBeadTypes(tools::Property &el) {
  for (auto &it : el) {
    if (it.name() == "rename") {
      string name = it.getAttribute<string>("name");
      string newname = it.getAttribute<string>("newname");
      if (name == "" || newname == "") {
        throw runtime_error("invalid rename tag, name or newname are empty.");
      }
      top_->RenameBeadType(name, newname);
    } else if (it.name() == "mass") {
      string name = it.getAttribute<string>("name");
      double value = it.getAttribute<double>("value");
      top_->SetBeadTypeMass(name, value);
    } else {
      throw std::runtime_error("Wrong element under beadtypes: " + it.name());
    }
  }
}

void XMLTopologyReader::ParseBonded(tools::Property &el) {
  for (auto &it : el) {
    if (it.name() == "bond") {
      ParseBond(it);
    } else if (it.name() == "angle") {
      ParseAngle(it);
    } else if (it.name() == "dihedral") {
      ParseDihedral(it);
    } else {
      throw std::runtime_error("Wrong element under bonded: " + it.name());
    }
  }
}

void XMLTopologyReader::ParseBond(tools::Property &p) {
  string name = p.get("name").as<string>();
  string beads = p.get("beads").as<string>();
  vector<string> bead_list = tools::Tokenizer(beads, " \n\t").ToVector();
  if (bead_list.size() % 2 == 1) {
    throw runtime_error("Wrong number of beads in bond: " + name);
  }
  Interaction *ic = nullptr;
  typedef pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  Index b_index = 0;
  for (vector<string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    if (b1.molname == b2.molname) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = molecules_.equal_range(b1.molname);
      for (MoleculesMap::iterator itm = mRange.first; itm != mRange.second;
           ++itm) {
        XMLMolecule &xmlMolecule = *itm->second;
        XMLBead &xmlBead1 = *xmlMolecule.name2beads[b1.atname];
        XMLBead &xmlBead2 = *xmlMolecule.name2beads[b2.atname];
        ic = new IBond(xmlBead1.pid, xmlBead2.pid);
        ic->setGroup(name);
        ic->setIndex(b_index);
        ic->setMolecule(xmlMolecule.pid);
        xmlMolecule.mi->AddInteraction(ic);
        top_->AddBondedInteraction(ic);
        b_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}

void XMLTopologyReader::ParseAngle(tools::Property &p) {
  string name = p.get("name").as<string>();
  string beads = p.get("beads").as<string>();
  vector<string> bead_list = tools::Tokenizer(beads, " \n\t").ToVector();
  if (bead_list.size() % 3 == 1) {
    throw runtime_error("Wrong number of beads in angle: " + name);
  }
  Interaction *ic = nullptr;
  typedef pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  Index b_index = 0;
  for (vector<string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    BondBead b3(*(it++));
    if ((b1.molname == b2.molname) && (b2.molname == b3.molname)) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = molecules_.equal_range(b1.molname);
      for (MoleculesMap::iterator itm = mRange.first; itm != mRange.second;
           ++itm) {
        XMLMolecule &xmlMolecule = *itm->second;
        XMLBead &xmlBead1 = *xmlMolecule.name2beads[b1.atname];
        XMLBead &xmlBead2 = *xmlMolecule.name2beads[b2.atname];
        XMLBead &xmlBead3 = *xmlMolecule.name2beads[b3.atname];
        ic = new IAngle(xmlBead1.pid, xmlBead2.pid, xmlBead3.pid);
        ic->setGroup(name);
        ic->setIndex(b_index);
        ic->setMolecule(xmlMolecule.pid);
        xmlMolecule.mi->AddInteraction(ic);
        top_->AddBondedInteraction(ic);
        b_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}
void XMLTopologyReader::ParseDihedral(tools::Property &p) {
  string name = p.get("name").as<string>();
  string beads = p.get("beads").as<string>();
  vector<string> bead_list = tools::Tokenizer(beads, " \n\t").ToVector();
  if (bead_list.size() % 4 == 1) {
    throw runtime_error("Wrong number of beads in dihedral: " + name);
  }
  Interaction *ic = nullptr;
  typedef pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  Index b_index = 0;
  for (vector<string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    BondBead b3(*(it++));
    BondBead b4(*(it++));
    if ((b1.molname == b2.molname) && (b3.molname == b4.molname) &&
        (b1.molname == b4.molname)) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = molecules_.equal_range(b1.molname);
      for (MoleculesMap::iterator itm = mRange.first; itm != mRange.second;
           ++itm) {
        XMLMolecule &xmlMolecule = *itm->second;
        XMLBead &xmlBead1 = *xmlMolecule.name2beads[b1.atname];
        XMLBead &xmlBead2 = *xmlMolecule.name2beads[b2.atname];
        XMLBead &xmlBead3 = *xmlMolecule.name2beads[b3.atname];
        XMLBead &xmlBead4 = *xmlMolecule.name2beads[b4.atname];
        ic = new IDihedral(xmlBead1.pid, xmlBead2.pid, xmlBead3.pid,
                           xmlBead4.pid);
        ic->setGroup(name);
        ic->setIndex(b_index);
        ic->setMolecule(xmlMolecule.pid);
        xmlMolecule.mi->AddInteraction(ic);
        top_->AddBondedInteraction(ic);
        b_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}

XMLTopologyReader::~XMLTopologyReader() {
  // Clean  molecules_ map
  for (auto &molecule_ : molecules_) {
    XMLMolecule *xmlMolecule = molecule_.second;
    for (auto &bead : xmlMolecule->beads) {
      delete bead;
    }
    xmlMolecule->beads.clear();
    delete xmlMolecule;
  }
}

}  // namespace csg
}  // namespace votca
