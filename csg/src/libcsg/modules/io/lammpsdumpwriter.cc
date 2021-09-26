/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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
#include <cstdio>
#include <string>
#include <votca/tools/constants.h>

namespace votca {
namespace csg {

using namespace std;
using namespace votca::tools;

void LAMMPSDumpWriter::Open(std::string file, bool bAppend) {
  out_ = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void LAMMPSDumpWriter::Close() { fclose(out_); }

void LAMMPSDumpWriter::Write(Topology *conf) {
  Topology *top = conf;
  Eigen::Matrix3d box = conf->getBox();
  fprintf(out_, "ITEM: TIMESTEP\n%ld\n", top->getStep());
  fprintf(out_, "ITEM: NUMBER OF ATOMS\n%li\n", (Index)top->Beads().size());
  fprintf(out_, "ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(out_, "0 %f\n0 %f\n0 %f\n", box(0, 0) * conv::nm2ang,
          box(1, 1) * conv::nm2ang, box(2, 2) * conv::nm2ang);

  fprintf(out_, "ITEM: ATOMS id type x y z");
  bool v = top->HasVel();
  if (v) {
    fprintf(out_, " vx vy vz");
  }
  bool f = top->HasForce();
  if (f) {
    fprintf(out_, " fx fy fz");
  }
  fprintf(out_, "\n");

  for (const Bead &bead : conf->Beads()) {
    Index type_id = conf->getBeadTypeId(bead.getType());

    fprintf(out_, "%ld %li", bead.getId() + 1, type_id);
    fprintf(out_, " %f %f %f", bead.getPos().x() * conv::nm2ang,
            bead.getPos().y() * conv::nm2ang, bead.getPos().z() * conv::nm2ang);
    if (v) {
      fprintf(out_, " %f %f %f", bead.getVel().x() * conv::nm2ang,
              bead.getVel().y() * conv::nm2ang,
              bead.getVel().z() * conv::nm2ang);
    }
    if (f) {
      fprintf(out_, " %f %f %f", bead.getF().x() * conv::kj2kcal / conv::nm2ang,
              bead.getF().y() * conv::kj2kcal / conv::nm2ang,
              bead.getF().z() * conv::kj2kcal / conv::nm2ang);
    }
    fprintf(out_, "\n");
  }
  fflush(out_);
}

}  // namespace csg
}  // namespace votca
