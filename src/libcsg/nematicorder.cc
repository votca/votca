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

#include <votca/csg/nematicorder.h>
#include <votca/csg/topology.h>
#include <votca/tools/matrix.h>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace csg {

using namespace std;

void NematicOrder::Process(Topology &top, const string &filter) {
  _mu.ZeroMatrix();
  _mv.ZeroMatrix();
  _mw.ZeroMatrix();
  int N = 0;
  bool bU, bV, bW;
  bU = bV = bW = false;

  for (BeadContainer::iterator iter = top.Beads().begin();
       iter != top.Beads().end(); ++iter) {

    Bead *bead = *iter;

    if (!wildcmp(filter.c_str(), bead->getName().c_str())) continue;

    if (bead->getSymmetry() == 1) continue;

    if (bead->HasU()) {
      _mu += bead->getU() | bead->getU();
      _mu[0][0] -= 1. / 3.;
      _mu[1][1] -= 1. / 3.;
      _mu[2][2] -= 1. / 3.;
      bU = true;
    }

    if (bead->HasV()) {
      _mv += bead->getV() | bead->getV();
      _mv[0][0] -= 1. / 3.;
      _mv[1][1] -= 1. / 3.;
      _mv[2][2] -= 1. / 3.;
      bV = true;
    }

    if (bead->HasW()) {
      _mw += bead->getW() | bead->getW();
      _mw[0][0] -= 1. / 3.;
      _mw[1][1] -= 1. / 3.;
      _mw[2][2] -= 1. / 3.;
      bW = true;
    }
    N++;
  }

  double f = 1. / (double)N * 3. / 2.;
  _mu = f * _mu;
  _mv = f * _mv;
  _mw = f * _mw;

  if (bU) _mu.SolveEigensystem(_nemat_u);
  if (bV) _mv.SolveEigensystem(_nemat_v);
  if (bW) _mw.SolveEigensystem(_nemat_w);
}

}  // namespace csg
}  // namespace votca
