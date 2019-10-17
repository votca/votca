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

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <votca/csg/topology.h>
#include <votca/tools/getline.h>
#include <votca/tools/tokenizer.h>

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include "dlpolytopologyreader.h"

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

string DLPOLYTopologyReader::_NextKeyline(ifstream &fs, const char *wspace)
// function to find and read the next line starting with a keyword/directive
// (skipping comments starting with "#" or ";") NOTE: the line is returned
// case-preserved, not to alter molecule/atom names (consider if --no-map is
// used)
{
  string line;
  size_t i_nws = 0;

  do {
    getline(fs, line);

    if (fs.eof()) {
      throw std::runtime_error("Error: unexpected end of dlpoly file '" +
                               _fname + "'");
    }

    i_nws = line.find_first_not_of(wspace);
  } while (line.substr(i_nws, 1) == "#" || line.substr(i_nws, 1) == ";");

  return line.substr(i_nws, line.size() - i_nws);
}

string DLPOLYTopologyReader::_NextKeyInt(ifstream &fs, const char *wspace,
                                         const string &word, int &ival)
// function to read the next line containing only a given keyword and an integer
// value after it (only skipping comments!) NOTE: this function must only be
// called when the next directive line has to contain the given keyword and an
// integer value
{
  stringstream sl(_NextKeyline(fs, wspace));
  string line, sval;

  sl >> line;  // allow user not to bother about the case
  boost::to_upper(line);

  if (line.substr(0, word.size()) != word) {
    throw std::runtime_error("Error: unexpected line from dlpoly file '" +
                             _fname + "', expected '" + word + "' but got '" +
                             line + "'");
  }

  sl >> sval;

  size_t i_num = sval.find_first_of(
      "0123456789");  // assume integer number straight after the only keyword

  if (i_num > 0) {
    throw std::runtime_error("Error: missing integer number in directive '" +
                             line + "' in topology file '" + _fname + "'");
  }

  ival = boost::lexical_cast<int>(sval);

  return sl.str();
}

bool DLPOLYTopologyReader::_isKeyInt(const string &line, const char *wspace,
                                     const string &word, int &ival)
// function to check if the given (last read) directive line starts with a given
// keyword and has an integer value at the end NOTE: comments are allowed beyond
// the integer (anything after the first integer is ignored)
{
  // split directives consisting of a few words: the sought keyword must be the
  // first one!
  Tokenizer tok(line, wspace);
  vector<string> fields;
  tok.ToVector(fields);

  ival = 0;

  if (fields.size() < 2) {
    return false;
  }

  boost::to_upper(fields[0]);

  if (fields[0].substr(0, word.size()) != word) {
    throw std::runtime_error("Error: unexpected directive from dlpoly file '" +
                             _fname + "', expected keyword '" + word +
                             "' but got '" + fields[0] + "'");
  }

  size_t i_num = string::npos;

  unsigned int i = 1;
  do {  // find integer number in the field with the lowest index (closest to
        // the keyword)
    i_num = fields[i++].find_first_of("0123456789");
  } while (i_num > 0 && i < fields.size());

  if (i_num > 0) {
    return false;
  }

  ival = boost::lexical_cast<int>(fields[i - 1]);

  return true;
}

bool DLPOLYTopologyReader::ReadTopology(string file, Topology &top) {
  const char *WhiteSpace = " \t";

  int matoms = 0;
  int natoms = 0;

  std::ifstream fl;
  boost::filesystem::path filepath(file.c_str());

  string line;

  // TODO: fix residue naming / assignment - DL_POLY has no means to recognise
  // residues!
  Residue *res = top.CreateResidue("no");

  if (boost::filesystem::basename(filepath).size() == 0) {
    if (filepath.parent_path().string().size() == 0) {
      _fname = "FIELD";  // DL_POLY uses fixed file names in current/working
                         // directory
    } else {
      _fname = filepath.parent_path().string() + "/FIELD";
    }
  } else {
    _fname = file;
  }

  fl.open(_fname.c_str());

  if (!fl.is_open()) {
    throw std::runtime_error("Error on opening dlpoly file '" + _fname + "'");
  } else {

    line = _NextKeyline(fl, WhiteSpace);  // read title line and skip it
    line = _NextKeyline(fl, WhiteSpace);  // read next directive line
    boost::to_upper(line);

    if (line.substr(0, 4) == "UNIT") {      // skip 'unit' line
      line = _NextKeyline(fl, WhiteSpace);  // read next directive line
      boost::to_upper(line);
    }

    if (line.substr(0, 4) == "NEUT") {  // skip 'neutral groups' line (DL_POLY
                                        // Classic FIELD format)
      line = _NextKeyline(fl, WhiteSpace);  // look for next directive line
      boost::to_upper(line);
    }

    int nmol_types;

    if (!_isKeyInt(line, WhiteSpace, "MOLEC", nmol_types)) {
      throw std::runtime_error("Error: missing integer number in directive '" +
                               line + "' in topology file '" + _fname + "'");
    }

#ifdef DEBUG
    cout << "Read from dlpoly file '" << _fname << "' : '" << line << "' - "
         << nmol_types << endl;
#endif

    string mol_name;

    for (int nmol_type = 0; nmol_type < nmol_types; nmol_type++) {

      mol_name = _NextKeyline(fl, WhiteSpace);
      Molecule *mi = top.CreateMolecule(mol_name);

      int nreplica = 1;
      line = _NextKeyInt(fl, WhiteSpace, "NUMMOL", nreplica);

#ifdef DEBUG
      cout << "Read from dlpoly file '" << _fname << "' : '" << mol_name
           << "' - '" << line << "' - " << nreplica << endl;
#endif

      line = _NextKeyInt(fl, WhiteSpace, "ATOMS", natoms);

#ifdef DEBUG
      cout << "Read from dlpoly file '" << _fname << "' : '" << line << "' - "
           << natoms << endl;
#endif

      // read molecule
      int id_map[natoms];
      for (int i = 0; i < natoms;) {  // i is altered in repeater loop
        stringstream sl(_NextKeyline(fl, WhiteSpace));

#ifdef DEBUG
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

        Tokenizer tok(line, WhiteSpace);
        vector<string> fields;
        tok.ToVector(fields);

#ifdef DEBUG
        cout << "Rest atom specs from dlpoly topology : '" << line << "'"
             << endl;
#endif

        int repeater = 1;
        if (fields.size() > 1) {
          repeater = boost::lexical_cast<int>(fields[0]);
        }

        for (int j = 0; j < repeater; j++) {

          string beadname = beadtype + "#" + boost::lexical_cast<string>(i + 1);
          Bead *bead =
              top.CreateBead(1, beadname, beadtype, res->getId(), mass, charge);

          stringstream nm;
          nm << bead->getResnr() + 1 << ":"
             << top.getResidue(bead->getResnr())->getName() << ":"
             << bead->getName();
          mi->AddBead(bead, nm.str());
          id_map[i] = bead->getId();
          i++;
#ifdef DEBUG
          cout << "Atom identification in maps '" << nm.str() << "'" << endl;
#endif
        }
        matoms += repeater;
      }

      while (line != "FINISH") {

        stringstream nl(_NextKeyline(fl, WhiteSpace));
        nl >> line;

#ifdef DEBUG
        cout << "Read unit type# from dlpoly topology : '" << nl.str() << "'"
             << endl;
#endif

        boost::to_upper(line);
        line = line.substr(0, 6);
        if ((line == "BONDS") || (line == "ANGLES") || (line == "DIHEDR")) {
          string type = line;
          int count;
          nl >> count;
          for (int i = 0; i < count; i++) {

            stringstream sl(_NextKeyline(fl, WhiteSpace));
#ifdef DEBUG
            cout << "Read unit specs from dlpoly topology : '" << sl.str()
                 << "'" << endl;
#endif
            sl >> line;  // internal dlpoly bond/angle/dihedral function types
                         // are merely skipped (ignored)
            int ids[4];
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

#ifdef DEBUG
      cout << "Read from dlpoly file '" << _fname << "' : '" << line
           << "' - done with '" << mol_name << "'" << endl;
#endif

      // replicate molecule
      for (int replica = 1; replica < nreplica; replica++) {
        Molecule *mi_replica = top.CreateMolecule(mol_name);
        for (int i = 0; i < mi->BeadCount(); i++) {
          Bead *bead = mi->getBead(i);
          string type = bead->getType();
          Bead *bead_replica =
              top.CreateBead(1, bead->getName(), type, res->getId(),
                             bead->getMass(), bead->getQ());
          mi_replica->AddBead(bead_replica, bead->getName());
        }
        matoms += mi->BeadCount();
        InteractionContainer ics = mi->Interactions();
        for (vector<Interaction *>::iterator ic = ics.begin(); ic != ics.end();
             ++ic) {
          Interaction *ic_replica = nullptr;
          int offset =
              mi_replica->getBead(0)->getId() - mi->getBead(0)->getId();
          if ((*ic)->BeadCount() == 2) {
            ic_replica = new IBond((*ic)->getBeadId(0) + offset,
                                   (*ic)->getBeadId(1) + offset);
          } else if ((*ic)->BeadCount() == 3) {
            ic_replica = new IAngle((*ic)->getBeadId(0) + offset,
                                    (*ic)->getBeadId(1) + offset,
                                    (*ic)->getBeadId(2) + offset);
          } else if ((*ic)->BeadCount() == 4) {
            ic_replica = new IDihedral(
                (*ic)->getBeadId(0) + offset, (*ic)->getBeadId(1) + offset,
                (*ic)->getBeadId(2) + offset, (*ic)->getBeadId(3) + offset);
          } else {
            throw std::runtime_error("Error: BeadCount not equal 2, 3 or 4");
          }
          ic_replica->setGroup((*ic)->getGroup());
          ic_replica->setIndex((*ic)->getIndex());
          ic_replica->setMolecule(mi_replica->getId());
          top.AddBondedInteraction(ic_replica);
          mi_replica->AddInteraction(ic_replica);
        }
      }
    }
    top.RebuildExclusions();
  }

#ifdef DEBUG
  getline(fl, line);  // is "close" found?
  if (line == "close") {
    cout << "Read from dlpoly file '" << _fname << "' : '" << line
         << "' - done with topology" << endl;
  } else {
    cout << "Read from dlpoly file '" << _fname
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
