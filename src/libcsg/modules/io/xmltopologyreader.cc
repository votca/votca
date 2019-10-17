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

#include "xmltopologyreader.h"
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>

namespace votca {
namespace csg {

using namespace std;

bool XMLTopologyReader::ReadTopology(string filename, Topology &top) {
  _top = &top;

  tools::Property options;
  options.LoadFromXML(filename);
  ParseRoot(options.get("topology"));

  _top->RebuildExclusions();
  return true;
}

void XMLTopologyReader::ReadTopolFile(string file) {
  TopologyReader *reader;
  reader = TopReaderFactory().Create(file);
  if (!reader) {
    throw runtime_error(file + ": unknown topology format");
  }

  reader->ReadTopology(file, *_top);

  delete reader;
  // Clean XML molecules and beads.
}

void XMLTopologyReader::ParseRoot(tools::Property &property) {
  _has_base_topology = false;
  if (property.hasAttribute("base")) {
    ReadTopolFile(property.getAttribute<string>("base"));
    _has_base_topology = true;
  }

  // Iterate over keys at first level.
  for (auto &it : property) {
    if (it.name() == "h5md_particle_group") {
      _top->setParticleGroup(it.getAttribute<string>("name"));
    } else if (it.name() == "molecules") {
      _mol_index = 1;
      _bead_index = 0;
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
  _top->setBox(m);
}

void XMLTopologyReader::ParseMolecules(tools::Property &p) {
  for (auto &it : p) {
    if (it.name() == "clear") {
      _top->ClearMoleculeList();
    } else if (it.name() == "rename") {
      string molname = it.getAttribute<string>("name");
      string range = it.getAttribute<string>("range");
      _top->RenameMolecules(range, molname);
    } else if (it.name() == "define" || it.name() == "molecule") {
      string molname = it.getAttribute<string>("name");
      int first = 0;
      if (it.name() == "define") {
        first = it.getAttribute<int>("first");
      }
      int nbeads = it.getAttribute<int>("nbeads");
      int nmols = it.getAttribute<int>("nmols");
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
        _top->CreateMoleculesByRange(molname, first, nbeads, nmols);
      } else {
        if (_has_base_topology) {
          throw std::runtime_error(
              "The defined list of beads only works for pure xml topology, "
              "without 'base' attribute.");
        }
        ParseMolecule(it, molname, nbeads, nmols);
      }
    }
  }
}

void XMLTopologyReader::ParseMolecule(tools::Property &p, string molname,
                                      int nbeads, int nmols) {
  vector<XMLBead *> xmlBeads;
  vector<int> xmlResidues;
  for (auto &it : p) {
    if (it.name() == "bead") {
      string atname = it.getAttribute<string>("name");
      string attype = it.getAttribute<string>("type");
      double atmass, atq;
      int resid;
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
        resid = it.getAttribute<int>("resid");
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
  int resnr = _top->ResidueCount();
  if (!xmlResidues.empty()) {
    if (xmlResidues.front() != resnr + 1 && xmlResidues.front() != -1) {
      throw std::runtime_error(
          "Residue count for beads in topology.molecules.molecule has to be "
          "greater than the number of residues already in the topology");
    }
  }
  for (int mn = 0; mn < nmols; mn++) {
    Molecule *mi = _top->CreateMolecule(molname);
    XMLMolecule *xmlMolecule = new XMLMolecule(molname, nmols);
    xmlMolecule->pid = mi->getId();
    xmlMolecule->mi = mi;
    _molecules.insert(make_pair(molname, xmlMolecule));
    vector<int>::iterator resit = xmlResidues.begin();
    for (vector<XMLBead *>::iterator itb = xmlBeads.begin();
         itb != xmlBeads.end(); ++itb, ++resit) {
      stringstream bname;
      XMLBead &b = **itb;
      if (*resit != -1) {
        if (_top->ResidueCount() < *resit) {
          resnr = *resit - 1;
          _top->CreateResidue(molname, resnr);
        }
      } else {
        _top->CreateResidue(molname, resnr);
      }

      if (!_top->BeadTypeExist(b.type)) {
        _top->RegisterBeadType(b.type);
      }
      Bead *bead = _top->CreateBead(1, b.name, b.type, resnr, b.mass, b.q);
      bname << _mol_index << ":" << molname << ":" << b.name;
      mi->AddBead(bead, bname.str());

      // Data for bonded terms.
      XMLBead *b_rep = new XMLBead(b);
      b_rep->pid = _bead_index;
      if (xmlMolecule->name2beads.count(b.name) != 0) {
        throw std::runtime_error("Atom " + b.name + " in molecule " + molname +
                                 " already exists.");
      }
      xmlMolecule->name2beads.insert(make_pair(b.name, b_rep));
      _bead_index++;
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
      _top->RenameBeadType(name, newname);
    } else if (it.name() == "mass") {
      string name = it.getAttribute<string>("name");
      double value = it.getAttribute<double>("value");
      _top->SetBeadTypeMass(name, value);
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
  tools::Tokenizer tok(beads, " \n\t");
  vector<string> bead_list = tok.ToVector();
  if (bead_list.size() % 2 == 1) {
    throw runtime_error("Wrong number of beads in bond: " + name);
  }
  Interaction *ic = nullptr;
  typedef pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  int b_index = 0;
  for (vector<string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    if (b1.molname == b2.molname) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = _molecules.equal_range(b1.molname);
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
        _top->AddBondedInteraction(ic);
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
  tools::Tokenizer tok(beads, " \n\t");
  vector<string> bead_list = tok.ToVector();
  if (bead_list.size() % 3 == 1) {
    throw runtime_error("Wrong number of beads in angle: " + name);
  }
  Interaction *ic = nullptr;
  typedef pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  int b_index = 0;
  for (vector<string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    BondBead b3(*(it++));
    if ((b1.molname == b2.molname) && (b2.molname == b3.molname)) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = _molecules.equal_range(b1.molname);
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
        _top->AddBondedInteraction(ic);
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
  tools::Tokenizer tok(beads, " \n\t");
  vector<string> bead_list = tok.ToVector();
  if (bead_list.size() % 4 == 1) {
    throw runtime_error("Wrong number of beads in dihedral: " + name);
  }
  Interaction *ic = nullptr;
  typedef pair<MoleculesMap::iterator, MoleculesMap::iterator> MRange;
  int b_index = 0;
  for (vector<string>::iterator it = bead_list.begin();
       it != bead_list.end();) {
    BondBead b1(*(it++));
    BondBead b2(*(it++));
    BondBead b3(*(it++));
    BondBead b4(*(it++));
    if ((b1.molname == b2.molname) && (b3.molname == b4.molname) &&
        (b1.molname == b4.molname)) {
      // Iterates over molecules and gets atom pids.
      MRange mRange = _molecules.equal_range(b1.molname);
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
        _top->AddBondedInteraction(ic);
        b_index++;
      }
    } else {
      throw std::runtime_error(
          "Beads from different molecules, not supported!");
    }
  }
}

XMLTopologyReader::~XMLTopologyReader() {
  // Clean _molecules map
  for (auto &_molecule : _molecules) {
    XMLMolecule *xmlMolecule = _molecule.second;
    for (auto &bead : xmlMolecule->beads) {
      delete bead;
    }
    xmlMolecule->beads.clear();
    delete xmlMolecule;
  }
}

}  // namespace csg
}  // namespace votca
