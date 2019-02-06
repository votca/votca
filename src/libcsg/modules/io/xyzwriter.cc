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

#include "xyzwriter.h"
#include <stdio.h>
#include <string>

namespace votca {
namespace csg {

using namespace std;

void XYZWriter::Open(string file, bool bAppend) {
  _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void XYZWriter::Close() { fclose(_out); }

void XYZWriter::Write(Topology *conf) {
  Topology *top = conf;
  fprintf(_out, "%d\n", (int)top->Beads().size());
  fprintf(_out, "frame: %d time: %f\n", top->getStep() + 1, top->getTime());

  for (BeadContainer::iterator iter = conf->Beads().begin();
       iter != conf->Beads().end(); ++iter) {
    Bead *bi = *iter;
    vec r = bi->getPos();
    // truncate strings if necessary
    string atomname = bi->getName();
    if (atomname.size() > 3) {
      atomname = atomname.substr(0, 3);
    }
    while (atomname.size() < 3) atomname = " " + atomname;

    // nm -> Angs
    fprintf(_out, "%s%10.5f%10.5f%10.5f\n", atomname.c_str(), r.getX() * 10.0,
            r.getY() * 10.0, r.getZ() * 10.0);
  }
  fflush(_out);
}

}  // namespace csg
}  // namespace votca
