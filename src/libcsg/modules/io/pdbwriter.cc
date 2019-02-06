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
#include <stdio.h>
#include <string>

namespace votca {
namespace csg {

using namespace std;

void PDBWriter::Open(string file, bool bAppend) {
  _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void PDBWriter::Close() { fclose(_out); }

void PDBWriter::Write(Topology *conf) {
  Topology *top = conf;
  fprintf(_out, "MODEL     %4d\n", conf->getStep());
  for (BeadContainer::iterator iter = conf->Beads().begin();
       iter != conf->Beads().end(); ++iter) {
    Bead *bi = *iter;
    vec r = bi->getPos();
    // truncate strings if necessary
    string resname = "";
    if (top->getResidue(bi->getResnr()))
      resname = top->getResidue(bi->getResnr())->getName();
    string atomname = bi->getName();
    if (resname.size() > 3) {
      resname = resname.substr(0, 3);
    }
    if (atomname.size() > 4) {
      atomname = atomname.substr(0, 4);
    }

    fprintf(_out, "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f\n",
            (bi->getId() + 1) % 100000,           // atom serial number
            atomname.c_str(),                     // atom name
            resname.c_str(),                      // residue name
            " ",                                  // chain identifier 1 char
            bi->getResnr() + 1,                   // residue sequence number
            10 * r.x(), 10 * r.y(), 10 * r.z());  // nm -> Angs
    // we skip the charge

    if (bi->getSymmetry() >= 2) {
      vec ru = 0.1 * bi->getU() + r;

      fprintf(_out, "HETATM%5d %4s %3s %1s%4d    %8.3f%8.3f%8.4f\n",
              bi->getId() + 1,          // atom serial number
              bi->getName().c_str(),    // atom name
              "REU",                    // residue name
              " ",                      // chain identifier 1 char
              bi->getResnr() + 1,       // residue sequence number
              ru.x(), ru.y(), ru.z());  // we skip the charge
    }
    if (bi->getSymmetry() >= 3) {
      vec rv = 0.1 * bi->getV() + r;
      fprintf(_out, "HETATM%5d %4s %3s %1s%4d    %8.3f%8.3f%8.4f\n",
              bi->getId() + 1,          // atom serial number
              bi->getName().c_str(),    // atom name
              "REV",                    // residue name
              " ",                      // chain identifier 1 char
              bi->getResnr() + 1,       // residue sequence number
              rv.x(), rv.y(), rv.z());  // we skip the charge
    }
  }
  fprintf(_out, "ENDMDL\n");
  fflush(_out);
}

}  // namespace csg
}  // namespace votca
