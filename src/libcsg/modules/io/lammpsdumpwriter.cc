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

#include "lammpsdumpwriter.h"
#include <stdio.h>
#include <string>

namespace votca {
namespace csg {

using namespace std;

void LAMMPSDumpWriter::Open(std::string file, bool bAppend) {
  _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void LAMMPSDumpWriter::Close() { fclose(_out); }

void LAMMPSDumpWriter::Write(Topology *conf) {
  Topology *      top = conf;
  Eigen::Matrix3d box = conf->getBox();
  fprintf(_out, "ITEM: TIMESTEP\n%i\n", top->getStep());
  fprintf(_out, "ITEM: NUMBER OF ATOMS\n%i\n", (int)top->Beads().size());
  fprintf(_out, "ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(_out, "0 %f\n0 %f\n0 %f\n", box(0, 0), box(1, 1), box(2, 2));

  fprintf(_out, "ITEM: ATOMS id type x y z");
  bool v = top->HasVel();
  if (v) {
    fprintf(_out, " vx vy vz");
  }
  bool f = top->HasForce();
  if (f) {
    fprintf(_out, " fx fy fz");
  }
  fprintf(_out, "\n");

  for (BeadContainer::iterator iter = conf->Beads().begin();
       iter != conf->Beads().end(); ++iter) {
    Bead *bi = *iter;

    int type_id = conf->getBeadTypeId(bi->getType());

    fprintf(_out, "%i %i", bi->getId() + 1, type_id);
    fprintf(_out, " %f %f %f", bi->getPos().x(), bi->getPos().y(),
            bi->getPos().z());
    if (v) {
      fprintf(_out, " %f %f %f", bi->getVel().x(), bi->getVel().y(),
              bi->getVel().z());
    }
    if (f) {
      fprintf(_out, " %f %f %f", bi->getF().x(), bi->getF().y(),
              bi->getF().z());
    }
    fprintf(_out, "\n");
  }
  fflush(_out);
}

}  // namespace csg
}  // namespace votca
