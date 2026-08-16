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

#pragma once
#ifndef VOTCA_XTP_FRAGMENTSATURATOR_H
#define VOTCA_XTP_FRAGMENTSATURATOR_H

// Local VOTCA includes
#include "qmmolecule.h"

namespace votca {
namespace xtp {

/**
 * \brief Saturates a mapped QMMolecule's own, real, external (i.e.
 * covalently bonded at the MD level, but mapped out of this
 * fragment/segment) bonds with new H atoms.
 *
 * This is deliberately the isolated, purely geometric placement step
 * only -- given a direction (already recorded on each QMAtom by
 * SegmentMapper, and already rigid-body-transformed into this
 * fragment's own, mapped/idealized reference frame -- see
 * QMAtom::hasExternalBond()/getExternalBondDirection() and
 * Md2QmEngine::map()'s own header comment for how that direction was
 * originally computed), places a new H atom at a standard, typical
 * bond length along it. It does NOT attempt any relaxation of the
 * new H's own position at all -- confirmed directly, this session, by
 * an actual, real DFT-level test, that this initial, purely geometric
 * placement alone carries a real, non-trivial error (~0.08 Angstrom,
 * for a deliberately, realistically distorted MD segment) -- a
 * follow-up, constrained relaxation step (planned separately, not yet
 * implemented) is expected to be necessary for most real use, not
 * merely an optional refinement. This class's own, single
 * responsibility is only the initial placement itself.
 */
class FragmentSaturator {
 public:
  // A generic, "typical" C-H bond length (Angstrom) -- not specific to
  // any one hybridization/element (this class only ever adds H atoms,
  // never any other element, matching the "H-saturation" scope this
  // whole feature was designed around from the start). Deliberately
  // just a reasonable, generic default, not a carefully-chosen,
  // atom-specific value -- the caller can override this directly via
  // SaturateExternalBonds's own bond_length_angstrom parameter, and,
  // in any case, this initial placement is only ever a starting point
  // for a later, separate relaxation step, not the final, trusted
  // position.
  static constexpr double kDefaultCHBondLengthAngstrom = 1.09;

  // Returns a NEW QMMolecule (the input mol is never modified) with an
  // additional H atom appended for every atom in mol that has
  // hasExternalBond() == true, placed at bond_length_angstrom along
  // that atom's own, already-recorded, already-transformed external-
  // bond direction. Atoms without an external bond are copied
  // unchanged, at their own, original indices -- the new H atoms are
  // appended strictly after all of mol's own, original atoms, so
  // every original atom's own index is preserved exactly (matching
  // the same, direct requirement already confirmed, earlier this
  // session, to be necessary for fragment_A/fragment_B atom-index
  // definitions computed before this saturation step to remain valid
  // afterward).
  static QMMolecule SaturateExternalBonds(
      const QMMolecule& mol,
      double bond_length_angstrom = kDefaultCHBondLengthAngstrom);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_FRAGMENTSATURATOR_H
