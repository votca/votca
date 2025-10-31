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

// Standard includes
#include <cstdio>
#include <string>

// Local private VOTCA includes
#include "growriter.h"

namespace votca {
namespace csg {

using namespace std;

void GROWriter::Open(string file, bool bAppend) {
  out_ = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void GROWriter::Close() { fclose(out_); }

void GROWriter::Write(Topology *conf) {

  constexpr int MAX_FMT_STR_LEN = 100;

  char format[MAX_FMT_STR_LEN];
  Index i, l, vpr;
  Topology *top = conf;

  fprintf(out_, "%s\n", "what a nice title");
  fprintf(out_, "%5ld\n", top->BeadCount());

  bool v = top->HasVel();
  Index pr = 3;  // precision of writeout, given by the spec

  /* build format string for printing,
     something like "%8.3f" for x and "%8.4f" for v */

  l = pr + 5;
  vpr = pr + 1;
  if (v) {
    snprintf(format, MAX_FMT_STR_LEN,
             "%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf\n",
             l, pr, l, pr, l, pr, l, vpr, l, vpr, l, vpr);
  } else {
    snprintf(format, MAX_FMT_STR_LEN, "%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf\n", l, pr,
             l, pr, l, pr);
  }

  for (i = 0; i < top->BeadCount(); i++) {
    Index resnr = top->getBead(i)->getResnr();
    string resname = top->getResidue(resnr).getName();
    string atomname = top->getBead(i)->getName();

    fprintf(out_, "%5ld%-5.5s%5.5s%5ld", (resnr + 1) % 100000, resname.c_str(),
            atomname.c_str(), (i + 1) % 100000);
    /* next fprintf uses built format string */
    Eigen::Vector3d r = conf->getBead(i)->getPos();

    if (v) {
      Eigen::Vector3d vv = conf->getBead(i)->getVel();
      fprintf(out_, format, r.x(), r.y(), r.z(), vv.x(), vv.y(), vv.z());
    } else {
      fprintf(out_, format, r.x(), r.y(), r.z());
    }
  }

  // write the box
  Eigen::Matrix3d box = conf->getBox();

  if (pr < 5) {
    pr = 5;
  }
  l = pr + 5;

  Eigen::Matrix3d box_offdiag = box;
  box_offdiag.diagonal().array() = 0.0;

  if (box_offdiag.isApproxToConstant(0, 1e-9)) {
    snprintf(format, MAX_FMT_STR_LEN,
             "%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf"
             "%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf\n",
             l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr);
    fprintf(out_, format, box(0, 0), box(1, 1), box(2, 2), box(1, 0), box(2, 0),
            box(0, 1), box(2, 1), box(0, 2), box(1, 2));
  } else {
    snprintf(format, MAX_FMT_STR_LEN, "%%%ld.%ldf%%%ld.%ldf%%%ld.%ldf\n", l, pr,
             l, pr, l, pr);
    fprintf(out_, format, box(0, 0), box(1, 1), box(2, 2));
  }
  fflush(out_);
}

}  // namespace csg
}  // namespace votca
