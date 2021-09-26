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

// Local VOTCA includes
#include "votca/csg/nematicorder.h"

namespace votca {
namespace csg {

using namespace std;

void NematicOrder::Process(Topology &top, const string &filter) {
  mu_ = Eigen::Matrix3d::Zero();
  mv_ = Eigen::Matrix3d::Zero();
  mw_ = Eigen::Matrix3d::Zero();
  Index N = 0;
  bool bU, bV, bW;
  bU = bV = bW = false;

  for (const auto &bead : top.Beads()) {

    if (!tools::wildcmp(filter, bead.getName())) {
      continue;
    }

    if (bead.getSymmetry() == 1) {
      continue;
    }

    if (bead.HasU()) {
      mu_ += bead.getU() * bead.getU().transpose();
      mu_.diagonal().array() -= 1. / 3.;
      bU = true;
    }

    if (bead.HasV()) {
      mu_ += bead.getV() * bead.getV().transpose();
      mu_.diagonal().array() -= 1. / 3.;
      bV = true;
    }

    if (bead.HasW()) {
      mu_ += bead.getW() * bead.getW().transpose();
      mu_.diagonal().array() -= 1. / 3.;
      bW = true;
    }
    N++;
  }

  double f = 1. / (double)N * 3. / 2.;
  mu_ = f * mu_;
  mv_ = f * mv_;
  mw_ = f * mw_;

  if (bU) {
    nemat_u_.computeDirect(mu_);
  }
  if (bV) {
    nemat_v_.computeDirect(mv_);
  }
  if (bW) {
    nemat_w_.computeDirect(mw_);
  }
}

}  // namespace csg
}  // namespace votca
