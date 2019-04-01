/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/aomatrix.h>
#include <votca/xtp/lowdin.h>

#include "votca/xtp/qmatom.h"
namespace votca {
namespace xtp {

void Lowdin::EvaluateLowdin(std::vector<QMAtom*>& _atomlist,
                            const Eigen::MatrixXd& _dmat, AOBasis& basis,
                            bool _do_transition) {
  AOOverlap _overlap;
  // Fill overlap
  _overlap.Fill(basis);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_overlap.Matrix());

  Eigen::MatrixXd Smsqrt = es.operatorSqrt();
  Eigen::MatrixXd _prodmat = Smsqrt * _dmat * Smsqrt;
  std::vector<QMAtom*>::iterator atom;

  int id = 0;
  for (atom = _atomlist.begin(); atom < _atomlist.end(); ++atom) {
    double charge = 0.0;
    // get element type and determine its nuclear charge
    if (!_do_transition) {
      charge = (*atom)->getNuccharge();
    }

    // a little messy, have to use basis set to find out which entries in dmat
    // belong to each atom.

    int nooffunc = basis.getFuncOfAtom((*atom)->getAtomID());

    for (int _i = id; _i < id + nooffunc; _i++) {
      charge -= _prodmat(_i, _i);
    }
    (*atom)->setPartialcharge(charge);
    id += nooffunc;
  }
  return;
}

}  // namespace xtp
}  // namespace votca
