/*
 * Copyright 2009-2023 The VOTCA Development Team (http://www.votca.org)
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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/getline.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/csg/topology.h"

#ifndef HAVE_NO_CONFIG
#include <votca_csg_config.h>
#endif

// Local private VOTCA includes
#include "dlpolytopologyreader.h"

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

string DLPOLYTopologyReader::NextKeyline_(ifstream &fs, const char *wspace)
// function to find and read the next line starting with a keyword/directive
// (skipping comments starting with "#" or ";") NOTE: the line is returned
// case-preserved, not to alter molecule/atom names (consider if --no-map is
// used)
{
  string line;
  size_t i_nws = 0;

  do {
    tools::getline(fs, line);

    if (fs.eof()) {
      throw std::runtime_error("Error: unexpected end of dlpoly file '" +
                               fname_ + "'");
    }

    i_nws = line.find_first_not_of(wspace);
  } while (line.substr(i_nws, 1) == "#" || line.substr(i_nws, 1) == ";");

  return line.substr(i_nws, line.size() - i_nws);
}

string DLPOLYTopologyReader::NextKeyInt_(ifstream &fs, const char *wspace,
                                         const string &word, Index &ival)
// function to read the next line containing only a given keyword and an integer
// value after it (only skipping comments!) NOTE: this function must only be
// called when the next directive line has to contain the given keyword and an
// integer value
{
  stringstream sl(NextKeyline_(fs, wspace));
  string line, sval;

  sl >> line;  // allow user not to bother about the case
  boost::to_upper(line);

  if (line.substr(0, word.size()) != word) {
    throw std::runtime_error("Error: unexpected line from dlpoly file '" +
                             fname_ + "', expected '" + word + "' but got '" +
                             line + "'");
  }

  sl >> sval;

  size_t i_num = sval.find_first_of(
      "0123456789");  // assume integer number straight after the only keyword

  if (i_num > 0) {
    throw std::runtime_error("Error: missing integer number in directive '" +
                             line + "' in topology file '" + fname_ + "'");
  }

  ival = boost::lexical_cast<Index>(sval);

  return sl.str();
}

bool DLPOLYTopologyReader::isKeyInt_(const string &line, const char *wspace,
                                     const string &word, Index &ival)
// function to check if the given (last read) directive line starts with a given
// keyword and has an integer value at the end NOTE: comments are allowed beyond
// the integer (anything after the first integer is ignored)
{
  // split directives consisting of a few words: the sought keyword must be the
  // first one!

  vector<string> fields = Tokenizer(line, wspace).ToVector();

  ival = 0;

  if (fields.size() < 2) {
    return false;
  }

  boost::to_upper(fields[0]);

  if (fields[0].substr(0, word.size()) != word) {
    throw std::runtime_error("Error: unexpected directive from dlpoly file '" +
                             fname_ + "', expected keyword '" + word +
                             "' but got '" + fields[0] + "'");
  }

  size_t i_num = string::npos;

  Index i = 1;
  do {  // find integer number in the field with the lowest index (closest to
        // the keyword)
    i_num = fields[i++].find_first_of("0123456789");
  } while (i_num > 0 && i < Index(fields.size()));

  if (i_num > 0) {
    return false;
  }

  ival = boost::lexical_cast<Index>(fields[i - 1]);

  return true;
}

bool DLPOLYTopologyReader::ReadTopology(string file, Topology &top) {

  Index natoms = 0;

  std::filesystem::path filepath(file.c_str());

  string line;

  // TODO: fix residue naming / assignment - DL_POLY has no means to recognise
  // residues!
  const Residue &res = top.CreateResidue("no");

  if (!filepath.has_stem()) {
    if (filepath.parent_path().string().size() == 0) {
      fname_ = "FIELD";  // DL_POLY uses fixed file names in current/working
                         // directory
    } else {
      fname_ = filepath.parent_path().string() + "/FIELD";
    }
  } else {
    fname_ = file;
  }
  std::ifstream fl;
  fl.open(fname_);

  if (!fl.is_open()) {
    throw std::runtime_error("Error on opening dlpoly file '" + fname_ + "'");
  } else {
    const char *WhiteSpace = " \t";
    NextKeyline_(fl, WhiteSpace);         // read title line and skip it
    line = NextKeyline_(fl, WhiteSpace);  // read next directive line
    boost::to_upper(line);

    if (line.substr(0, 4) == "UNIT") {      // skip 'unit' line
      line = NextKeyline_(fl, WhiteSpace);  // read next directive line
      boost::to_upper(line);
    }

    if (line.substr(0, 4) == "NEUT") {  // skip 'neutral groups' line (DL_POLY
                                        // Classic FIELD format)
      line = NextKeyline_(fl, WhiteSpace);  // look for next directive line
      boost::to_upper(line);
    }

    Index nmol_types;

    if (!isKeyInt_(line, WhiteSpace, "MOLEC", nmol_types)) {
      throw std::runtime_error("Error: missing integer number in directive '" +
                               line + "' in topology file '" + fname_ + "'");
    }

#ifndef NDEBUG
    cout << "Read from dlpoly file '" << fname_ << "' : '" << line << "' - "
         << nmol_types << endl;
#endif

    string mol_name;

    for (Index nmol_type = 0; nmol_type < nmol_types; nmol_type++) {

      mol_name = NextKeyline_(fl, WhiteSpace);
      Molecule *mi = top.CreateMolecule(mol_name);

      Index nreplica = 1;
      line = NextKeyInt_(fl, WhiteSpace, "NUMMOL", nreplica);

#ifndef NDEBUG
      cout << "Read from dlpoly file '" << fname_ << "' : '" << mol_name
           << "' - '" << line << "' - " << nreplica << endl;
#endif

      line = NextKeyInt_(fl, WhiteSpace, "ATOMS", natoms);

#ifndef NDEBUG
      cout << "Read from dlpoly file '" << fname_ << "' : '" << line << "' - "
           << natoms << endl;
#endif
      std::vector<Index> id_map(natoms);
      for (Index i = 0; i < natoms;) {  // i is altered in repeater loop
        stringstream sl(NextKeyline_(fl, WhiteSpace));

#ifndef NDEBUG
        cout << "Read atom specs from dlpoly topology : '" << sl.str() << "'"
             << endl;
#endif
        string beadtype;
        sl >> beadtype;
        if (!top.BeadTypeExist(beadtype)) {
          top.RegisterBeadType(beadtype);
        }
        double mass;
        sl >> mass;
        double charge;
        sl >> charge;

        line = " ";
        sl >> line;  // rest of the atom line

        vector<string> fields = Tokenizer(line, WhiteSpace).ToVector();

#ifndef NDEBUG
        cout << "Rest atom specs from dlpoly topology : '" << line << "'"
             << endl;
#endif

        Index repeater = 1;
        if (fields.size() > 1) {
          repeater = boost::lexical_cast<Index>(fields[0]);
        }

        for (Index j = 0; j < repeater; j++) {

          string beadname = beadtype + "#" + boost::lexical_cast<string>(i + 1);
          Bead *bead = top.CreateBead(Bead::spherical, beadname, beadtype,
                                      res.getId(), mass, charge);

          stringstream nm;
          nm << bead->getResnr() + 1 << ":"
             << top.getResidue(bead->getResnr()).getName() << ":"
             << bead->getName();
          mi->AddBead(bead, nm.str());
          id_map[i] = bead->getId();
          i++;
#ifndef NDEBUG
          cout << "Atom identification in maps '" << nm.str() << "'" << endl;
#endif
        }
      }

      while (line != "FINISH") {

        stringstream nl(NextKeyline_(fl, WhiteSpace));
        nl >> line;

#ifndef NDEBUG
        cout << "Read unit type# from dlpoly topology : '" << nl.str() << "'"
             << endl;
#endif

        boost::to_upper(line);
        line = line.substr(0, 6);
        if ((line == "BONDS") || (line == "ANGLES") || (line == "DIHEDR")) {
          string type = line;
          Index count;
          nl >> count;
          for (Index i = 0; i < count; i++) {

            stringstream sl(NextKeyline_(fl, WhiteSpace));
#ifndef NDEBUG
            cout << "Read unit specs from dlpoly topology : '" << sl.str()
                 << "'" << endl;
#endif
            sl >> line;  // internal dlpoly bond/angle/dihedral function types
                         // are merely skipped (ignored)
            Index ids[4];
            Interaction *ic = nullptr;
            sl >> ids[0];
            sl >> ids[1];
            if (type == "BONDS") {
              ic = new IBond(id_map[ids[0] - 1],
                             id_map[ids[1] - 1]);  // -1 due to fortran vs c
            } else if (type == "ANGLES") {
              sl >> ids[2];
              ic = new IAngle(id_map[ids[0] - 1], id_map[ids[1] - 1],
                              id_map[ids[2] - 1]);  // -1 due to fortran vs c
            } else if (type.substr(0, 6) == "DIHEDR") {
              type = "DIHEDRALS";
              sl >> ids[2];
              sl >> ids[3];
              ic = new IDihedral(id_map[ids[0] - 1], id_map[ids[1] - 1],
                                 id_map[ids[2] - 1],
                                 id_map[ids[3] - 1]);  // -1 due to fortran vs c
            } else {
              throw std::runtime_error(
                  "Error: type should be BONDS, ANGLES or DIHEDRALS");
            }
            // could one use bond/angle/dihedral function types for 1:1 mapping?
            // (CG map overwrites ic->Group anyway)
            // ic->setGroup(line);
            ic->setGroup(type);
            ic->setIndex(i);
            ic->setMolecule(mi->getId());
            top.AddBondedInteraction(ic);
            mi->AddInteraction(ic);
          }
        }
      }

#ifndef NDEBUG
      cout << "Read from dlpoly file '" << fname_ << "' : '" << line
           << "' - done with '" << mol_name << "'" << endl;
#endif

      // replicate molecule
      for (Index replica = 1; replica < nreplica; replica++) {
        Molecule *mi_replica = top.CreateMolecule(mol_name);
        for (Index i = 0; i < mi->BeadCount(); i++) {
          Bead *bead = mi->getBead(i);
          string type = bead->getType();
          Bead *bead_replica =
              top.CreateBead(Bead::spherical, bead->getName(), type,
                             res.getId(), bead->getMass(), bead->getQ());
          mi_replica->AddBead(bead_replica, bead->getName());
        }
        InteractionContainer ics = mi->Interactions();
        for (auto &ic : ics) {
          Interaction *ic_replica = nullptr;
          Index offset =
              mi_replica->getBead(0)->getId() - mi->getBead(0)->getId();
          if (ic->BeadCount() == 2) {
            ic_replica =
                new IBond(ic->getBeadId(0) + offset, ic->getBeadId(1) + offset);
          } else if (ic->BeadCount() == 3) {
            ic_replica =
                new IAngle(ic->getBeadId(0) + offset, ic->getBeadId(1) + offset,
                           ic->getBeadId(2) + offset);
          } else if (ic->BeadCount() == 4) {
            ic_replica = new IDihedral(
                ic->getBeadId(0) + offset, ic->getBeadId(1) + offset,
                ic->getBeadId(2) + offset, ic->getBeadId(3) + offset);
          } else {
            throw std::runtime_error("Error: BeadCount not equal 2, 3 or 4");
          }
          ic_replica->setGroup(ic->getGroup());
          ic_replica->setIndex(ic->getIndex());
          ic_replica->setMolecule(mi_replica->getId());
          top.AddBondedInteraction(ic_replica);
          mi_replica->AddInteraction(ic_replica);
        }
      }
    }
    top.RebuildExclusions();
  }

#ifndef NDEBUG
  tools::getline(fl, line);  // is "close" found?
  if (line == "close") {
    cout << "Read from dlpoly file '" << fname_ << "' : '" << line
         << "' - done with topology" << endl;
  } else {
    cout << "Read from dlpoly file '" << fname_
         << "' : 'EOF' - done with topology (directive 'close' not read!)"
         << endl;
  }
#endif

  // we don't need the rest
  fl.close();

  return true;
}

}  // namespace csg
}  // namespace votca
