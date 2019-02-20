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

#include "growriter.h"
#include <stdio.h>
#include <string>

namespace votca {
namespace csg {

using namespace std;

void GROWriter::Open(string file, bool bAppend) {
  _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void GROWriter::Close() { fclose(_out); }

void GROWriter::Write(Topology *conf) {
  char format[100];
  int i, resnr, l, vpr;
  Topology *top = conf;

  fprintf(_out, "%s\n", "what a nice title");
  fprintf(_out, "%5d\n", top->BeadCount());

  bool v = top->HasVel();
  int pr = 3; // precision of writeout, given by the spec

  /* build format string for printing,
     something like "%8.3f" for x and "%8.4f" for v */
  /*if (pr<0)
    pr=0;
  if (pr>30)
    pr=30;*/
  l = pr + 5;
  vpr = pr + 1;
  if (v)
    sprintf(format, "%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n", l, pr,
            l, pr, l, pr, l, vpr, l, vpr, l, vpr);
  else
    sprintf(format, "%%%d.%df%%%d.%df%%%d.%df\n", l, pr, l, pr, l, pr);

  for (i = 0; i < top->BeadCount(); i++) {
    resnr = top->getBead(i)->getResnr();
    string resname = top->getResidue(resnr)->getName();
    string atomname = top->getBead(i)->getName();

    fprintf(_out, "%5d%-5.5s%5.5s%5d", (resnr + 1) % 100000, resname.c_str(),
            atomname.c_str(), (i + 1) % 100000);
    /* next fprintf uses built format string */
    Eigen::Vector3d r = conf->getBead(i)->getPos();

    if (v) {
      Eigen::Vector3d vv = conf->getBead(i)->getVel();
      fprintf(_out, format, r.x(), r.y(), r.z(), vv.x(), vv.y(), vv.z());
    } else {
      fprintf(_out, format, r.x(), r.y(), r.z());
    }
  }

  // write the boy
  Eigen::Matrix3d box = conf->getBox();

  if (pr < 5)
    pr = 5;
  l = pr + 5;

  if (box(0, 1) || box(0, 2) || box(1, 0) || box(1, 2) || box(2, 0) ||
      box(2, 1)) {
    sprintf(format,
            "%%%d.%df%%%d.%df%%%d.%df"
            "%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
            l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr);
    fprintf(_out, format, box(0, 0), box(1, 1), box(2, 2), box(1, 0), box(2, 0),
            box(0, 1), box(2, 1), box(0, 2), box(1, 2));
  } else {
    sprintf(format, "%%%d.%df%%%d.%df%%%d.%df\n", l, pr, l, pr, l, pr);
    fprintf(_out, format, box(0, 0), box(1, 1), box(2, 2));
  }
  fflush(_out);
}

} // namespace csg
} // namespace votca
