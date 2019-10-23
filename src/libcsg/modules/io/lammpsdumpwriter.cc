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
#include <votca/tools/constants.h>

namespace votca {
namespace csg {

using namespace std;
using namespace votca::tools;

void LAMMPSDumpWriter::Open(std::string file, bool bAppend) {
  _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void LAMMPSDumpWriter::Close() { fclose(_out); }

void LAMMPSDumpWriter::Write(Topology *conf) {
  Topology *top = conf;
  Eigen::Matrix3d box = conf->getBox();
  fprintf(_out, "ITEM: TIMESTEP\n%ld\n", top->getStep());
  fprintf(_out, "ITEM: NUMBER OF ATOMS\n%i\n", (int)top->Beads().size());
  fprintf(_out, "ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(_out, "0 %f\n0 %f\n0 %f\n", box(0, 0) * conv::nm2ang,
          box(1, 1) * conv::nm2ang, box(2, 2) * conv::nm2ang);

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

    long int type_id = conf->getBeadTypeId(bi->getType());

    fprintf(_out, "%ld %li", bi->getId() + 1, type_id);
    fprintf(_out, " %f %f %f", bi->getPos().x() * conv::nm2ang,
            bi->getPos().y() * conv::nm2ang, bi->getPos().z() * conv::nm2ang);
    if (v) {
      fprintf(_out, " %f %f %f", bi->getVel().x() * conv::nm2ang,
              bi->getVel().y() * conv::nm2ang, bi->getVel().z() * conv::nm2ang);
    }
    if (f) {
      fprintf(_out, " %f %f %f", bi->getF().x() * conv::kj2kcal / conv::nm2ang,
              bi->getF().y() * conv::kj2kcal / conv::nm2ang,
              bi->getF().z() * conv::kj2kcal / conv::nm2ang);
    }
    fprintf(_out, "\n");
  }
  fflush(_out);
}

}  // namespace csg
}  // namespace votca
