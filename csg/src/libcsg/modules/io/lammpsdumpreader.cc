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
#include <memory>
#include <vector>

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/getline.h>

// Local private VOTCA includes
#include "lammpsdumpreader.h"

namespace votca {
namespace csg {
using namespace boost;
using namespace std;

bool LAMMPSDumpReader::ReadTopology(string file, Topology &top) {
  topology_ = true;
  top.Cleanup();

  fl_.open(file);
  if (!fl_.is_open()) {
    throw std::ios_base::failure("Error on open topology file: " + file);
  }
  fname_ = file;

  NextFrame(top);

  fl_.close();

  return true;
}

bool LAMMPSDumpReader::Open(const string &file) {
  fl_.open(file);
  if (!fl_.is_open()) {
    throw std::ios_base::failure("Error on open trajectory file: " + file);
  }
  fname_ = file;
  return true;
}

void LAMMPSDumpReader::Close() { fl_.close(); }

bool LAMMPSDumpReader::FirstFrame(Topology &top) {
  topology_ = false;
  NextFrame(top);
  return true;
}

bool LAMMPSDumpReader::NextFrame(Topology &top) {
  string line;
  tools::getline(fl_, line);
  boost::algorithm::trim(line);
  while (!fl_.eof()) {
    if (line.substr(0, 5) != "ITEM:") {
      throw std::ios_base::failure("unexpected line in lammps file:\n" + line);
    }
    if (line.substr(6, 8) == "TIMESTEP") {
      ReadTimestep(top);
    } else if (line.substr(6, 15) == "NUMBER OF ATOMS") {
      ReadNumAtoms(top);
    } else if (line.substr(6, 10) == "BOX BOUNDS") {
      ReadBox(top);
    } else if (line.substr(6, 5) == "ATOMS") {
      ReadAtoms(top, line);
      break;
    }

    else {
      throw std::ios_base::failure("unknown item lammps file : " +
                                   line.substr(6));
    }
    tools::getline(fl_, line);
    boost::algorithm::trim(line);
  }
  if (topology_) {
    cout << "WARNING: topology created from .dump file, masses, charges, "
            "types, residue names are wrong!\n";
  }
  return !fl_.eof();
  ;
}

void LAMMPSDumpReader::ReadTimestep(Topology &top) {
  string s;
  tools::getline(fl_, s);
  boost::algorithm::trim(s);
  top.setStep(boost::lexical_cast<Index>(s));
  cout << "Reading frame, timestep " << top.getStep() << endl;
}

void LAMMPSDumpReader::ReadBox(Topology &top) {
  string s;

  Eigen::Matrix3d m = Eigen::Matrix3d::Zero();

  for (Index i = 0; i < 3; ++i) {
    tools::getline(fl_, s);
    boost::algorithm::trim(s);
    vector<double> v = tools::Tokenizer(s, " ").ToVector<double>();
    if (v.size() != 2) {
      throw std::ios_base::failure("invalid box format");
    }
    m(i, i) = v[1] - v[0];
  }
  top.setBox(m * tools::conv::ang2nm);
}

void LAMMPSDumpReader::ReadNumAtoms(Topology &top) {
  string s;
  tools::getline(fl_, s);
  boost::algorithm::trim(s);
  natoms_ = boost::lexical_cast<Index>(s);
  if (!topology_ && natoms_ != top.BeadCount()) {
    std::runtime_error("number of beads in topology and trajectory differ");
  }
}

void LAMMPSDumpReader::ReadAtoms(Topology &top, string itemline) {
  if (topology_) {
    top.CreateResidue("dum");
    if (!top.BeadTypeExist("no")) {
      top.RegisterBeadType("no");
    }
    for (Index i = 0; i < natoms_; ++i) {
      (void)top.CreateBead(Bead::spherical, "no", "no", 0, 0, 0);
    }
  }

  bool pos = false;
  bool force = false;
  bool vel = false;
  Index id = -1;

  vector<string> fields;

  {
    tools::Tokenizer tok(itemline.substr(12), " ");
    fields = tok.ToVector();
    Index j = 0;
    for (tools::Tokenizer::iterator i = tok.begin(); i != tok.end(); ++i, ++j) {
      if (*i == "x" || *i == "y" || *i == "z") {
        pos = true;
      } else if (*i == "xu" || *i == "yu" || *i == "zu") {
        pos = true;
      } else if (*i == "xs" || *i == "ys" || *i == "zs") {
        pos = true;
      } else if (*i == "vx" || *i == "vy" || *i == "vz") {
        vel = true;
      } else if (*i == "fx" || *i == "fy" || *i == "fz") {
        force = true;
      } else if (*i == "id") {
        id = j;
      }
    }
  }
  if (id < 0) {
    throw std::runtime_error(
        "error, id not found in any column of the atoms section");
  }

  for (Index i = 0; i < natoms_; ++i) {
    string s;
    tools::getline(fl_, s);
    boost::algorithm::trim(s);
    if (fl_.eof()) {
      throw std::runtime_error("Error: unexpected end of lammps file '" +
                               fname_ + "' only " +
                               boost::lexical_cast<string>(i) + " atoms of " +
                               boost::lexical_cast<string>(natoms_) + " read.");
    }

    tools::Tokenizer tok(s, " ");
    tools::Tokenizer::iterator itok = tok.begin();
    vector<string> fields2 = tok.ToVector();
    // internal numbering begins with 0
    Index atom_id = boost::lexical_cast<Index>(fields2[id]);
    if (atom_id > natoms_) {
      throw std::runtime_error(
          "Error: found atom with id " + boost::lexical_cast<string>(atom_id) +
          " but only " + boost::lexical_cast<string>(natoms_) +
          " atoms defined in header of file '" + fname_ + "'");
    }
    Bead *b = top.getBead(atom_id - 1);
    b->HasPos(pos);
    b->HasF(force);
    b->HasVel(vel);
    Eigen::Matrix3d m = top.getBox();

    for (size_t j = 0; itok != tok.end(); ++itok, ++j) {
      if (j == fields.size()) {
        throw std::runtime_error(
            "error, wrong number of columns in atoms section");
      } else if (fields[j] == "x") {
        b->Pos().x() = stod(*itok) * tools::conv::ang2nm;
      } else if (fields[j] == "y") {
        b->Pos().y() = stod(*itok) * tools::conv::ang2nm;
      } else if (fields[j] == "z") {
        b->Pos().z() = stod(*itok) * tools::conv::ang2nm;
      } else if (fields[j] == "xu") {
        b->Pos().x() = stod(*itok) * tools::conv::ang2nm;
      } else if (fields[j] == "yu") {
        b->Pos().y() = stod(*itok) * tools::conv::ang2nm;
      } else if (fields[j] == "zu") {
        b->Pos().z() = stod(*itok) * tools::conv::ang2nm;
      } else if (fields[j] == "xs") {
        b->Pos().x() = stod(*itok) * m(0, 0);  // box is already in nm
      } else if (fields[j] == "ys") {
        b->Pos().y() = stod(*itok) * m(1, 1);  // box is already in nm
      } else if (fields[j] == "zs") {
        b->Pos().z() = stod(*itok) * m(2, 2);  // box is already in nm
      } else if (fields[j] == "vx") {
        b->Vel().x() = stod(*itok) * tools::conv::ang2nm;
      } else if (fields[j] == "vy") {
        b->Vel().y() = stod(*itok) * tools::conv::ang2nm;
      } else if (fields[j] == "vz") {
        b->Vel().z() = stod(*itok) * tools::conv::ang2nm;
      } else if (fields[j] == "fx") {
        b->F().x() = stod(*itok) * tools::conv::kcal2kj / tools::conv::ang2nm;
      } else if (fields[j] == "fy") {
        b->F().y() = stod(*itok) * tools::conv::kcal2kj / tools::conv::ang2nm;
      } else if (fields[j] == "fz") {
        b->F().z() = stod(*itok) * tools::conv::kcal2kj / tools::conv::ang2nm;
      } else if ((fields[j] == "type") && topology_) {
        if (!top.BeadTypeExist(*itok)) {
          top.RegisterBeadType(*itok);
        }
        b->setType(*itok);
      }
    }
  }
}

}  // namespace csg
}  // namespace votca
