/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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
#include <string>

// Local private VOTCA includes
#include "gmxtrajectorywriter.h"

// Third party includes
#include <gromacs/fileio/trxio.h>
#include <gromacs/trajectory/trajectoryframe.h>

// this one is needed because of bool is defined in one of the headers included
// by gmx
#undef bool

namespace votca {
namespace csg {

void GMXTrajectoryWriter::Open(std::string file, bool) {
  file_ = open_trx((char *)file.c_str(), "w");
}

void GMXTrajectoryWriter::Close() { close_trx(file_); }

void GMXTrajectoryWriter::Write(Topology *conf) {
  Index N = conf->BeadCount();
  t_trxframe frame;
  rvec *x = new rvec[N];
  rvec *v = nullptr;
  rvec *f = nullptr;
  Eigen::Matrix3d box = conf->getBox();

  frame.natoms = (int)N;
  frame.bTime = true;
  frame.time = real(conf->getTime());
  frame.bStep = true;
  frame.x = x;
  frame.bLambda = false;
  frame.bAtoms = false;
  frame.bPrec = false;
  frame.bX = true;
  frame.bF = conf->HasForce();
  frame.bBox = true;
  frame.bV = conf->HasVel();

  for (Index i = 0; i < 3; i++) {
    for (Index j = 0; j < 3; j++) {
      frame.box[j][i] = real(box(i, j));
    }
  }

  for (Index i = 0; i < N; ++i) {
    Eigen::Vector3d pos = conf->getBead(i)->getPos();
    x[i][0] = real(pos.x());
    x[i][1] = real(pos.y());
    x[i][2] = real(pos.z());
  }

  if (frame.bV) {
    v = new rvec[N];
    for (Index i = 0; i < N; ++i) {
      frame.v = v;
      Eigen::Vector3d vel = conf->getBead(i)->getVel();
      v[i][0] = real(vel.x());
      v[i][1] = real(vel.y());
      v[i][2] = real(vel.z());
    }
  }
  if (frame.bF) {
    f = new rvec[N];
    for (Index i = 0; i < N; ++i) {
      frame.f = f;
      Eigen::Vector3d force = conf->getBead(i)->getF();
      f[i][0] = real(force.x());
      f[i][1] = real(force.y());
      f[i][2] = real(force.z());
    }
  }

  write_trxframe(file_, &frame, nullptr);

  delete[] x;
  if (frame.bV) {
    delete[] v;
  }
  if (frame.bF) {
    delete[] f;
  }
}

}  // namespace csg
}  // namespace votca
