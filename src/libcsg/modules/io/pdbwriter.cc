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

#include "pdbwriter.h"
#include <boost/format.hpp>
#include <stdio.h>
#include <string>

namespace votca {
namespace csg {

using namespace std;

void PDBWriter::Open(string file, bool bAppend) {
  if (bAppend) {
    _out.open(file, std::ios_base::app);
  } else {
    _out.open(file);
  }
}

void PDBWriter::Close() { _out.close(); }

void PDBWriter::Write(Topology *conf) {
  _out << boost::format("MODEL     %1$4d\n") % conf->getStep();
  boost::format atomfrmt(
      "ATOM  %1$5d %2$4s %3$3s %4$1s%5$4d    %6$8.3f%7$8.3f%8$8.3f\n");
  boost::format beadfrmt(
      "HETATM%1$5d %2$4s %3$3s %4$1s%5$4d    %6$8.3f%7$8.3f%8$8.3f\n");
  for (Bead *bi : conf->Beads()) {
    vec r = 10 * bi->getPos(); // nm -> Angs
    // truncate strings if necessary
    string resname = "";
    if (conf->getResidue(bi->getResnr()))
      resname = conf->getResidue(bi->getResnr())->getName();
    string atomname = bi->getName();
    if (resname.size() > 3) {
      resname = resname.substr(0, 3);
    }
    if (atomname.size() > 4) {
      atomname = atomname.substr(0, 4);
    }

    _out << atomfrmt % ((bi->getId() + 1) % 100000) // atom serial number
                % atomname % resname % " "          // chain identifier 1 char
                % (bi->getResnr() + 1)              // residue sequence number
                % r.x() % r.y() % r.z();
    // we skip the charge

    if (bi->getSymmetry() >= 2) {
      vec ru = 0.1 * bi->getU() + r;

      _out << beadfrmt % (bi->getId() + 1) % 100000 // atom serial number
                  % bi->getName()                   // atom name
                  % "REU"                           // residue name
                  % " "                             // chain identifier 1 char
                  % (bi->getResnr() + 1)            // residue sequence number
                  % ru.x() % ru.y() % ru.z();       // we skip the charge
    }
    if (bi->getSymmetry() >= 3) {
      vec rv = 0.1 * bi->getV() + r;
      _out << beadfrmt % (bi->getId() + 1) % 100000 // atom serial number
                  % bi->getName()                   // atom name
                  % "REV"                           // residue name
                  % " "                             // chain identifier 1 char
                  % (bi->getResnr() + 1)            // residue sequence number
                  % rv.x() % rv.y() % rv.z();       // we skip the charge
    }
  }
  _out << "ENDMDL\n";
  _out << std::flush;
}
}
}
