/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/fragmentsaturator.h"

namespace votca {
namespace xtp {

QMMolecule FragmentSaturator::SaturateExternalBonds(
    const QMMolecule& mol, double bond_length_angstrom) {
  QMMolecule result(mol.getType(), mol.getId());

  // All of mol's own, original atoms are copied first, entirely
  // unchanged, at their own, original indices -- new H atoms are only
  // ever appended after all of these, so every original index remains
  // valid, matching the same, direct requirement already confirmed,
  // earlier this session, for fragment_A/fragment_B atom-index
  // definitions computed before this saturation step.
  for (const QMAtom& atom : mol) {
    result.push_back(atom);
  }

  double bond_length_bohr = bond_length_angstrom * tools::conv::ang2bohr;
  Index new_index = mol.size();
  for (const QMAtom& atom : mol) {
    if (!atom.hasExternalBond()) {
      continue;
    }
    // getExternalBondDirection() is already normalized (confirmed
    // directly: QMAtom::setExternalBondDirection's own .normalize()
    // call, and a pure rotation, both preserve unit length exactly --
    // SegmentMapper's own transform never renormalizes this
    // separately, since it does not need to).
    Eigen::Vector3d new_pos =
        atom.getPos() + bond_length_bohr * atom.getExternalBondDirection();
    QMAtom new_h(new_index, "H", new_pos);
    result.push_back(new_h);
    new_index++;
  }

  return result;
}

}  // namespace xtp
}  // namespace votca
