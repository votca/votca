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
 * \brief Saturates the external bonds of a mapped QMMolecule (bonds
 * that exist at the MD level but cross the segment boundary) with H
 * atoms.
 *
 * SaturateExternalBonds places each H at a fixed bond length along the
 * external-bond direction stored on the QMAtom by SegmentMapper (already
 * transformed into the mapped frame; see
 * QMAtom::getExternalBondDirection() and Md2QmEngine::map()). This
 * placement alone has errors of about 0.08 Angstrom for distorted MD
 * geometries, so RelaxNewAtoms subsequently relaxes the added H atoms
 * with all other atoms fixed.
 */
class FragmentSaturator {
 public:
  // Generic C-H bond length [Angstrom], independent of hybridization.
  // Only the starting point for RelaxNewAtoms.
  static constexpr double kDefaultCHBondLengthAngstrom = 1.09;

  // new_atom_parent_ids has size mol.size() and is indexed by atom id:
  // -1 for original atoms, and for each added H the id of the original
  // atom it saturates. Callers use it to assign added H atoms to the
  // fragment of their parent (e.g. IPodCoupling::EvalJob); the saturator
  // itself has no notion of fragments.
  struct SaturationResult {
    QMMolecule mol;
    std::vector<Index> new_atom_parent_ids;
  };

  // Returns a copy of mol with one H atom appended for every atom with
  // hasExternalBond(), placed at bond_length_angstrom along its external
  // bond direction. Original atoms keep their indices, so atom indices
  // computed before saturation remain valid.
  static SaturationResult SaturateExternalBonds(
      const QMMolecule& mol,
      double bond_length_angstrom = kDefaultCHBondLengthAngstrom);
  // Returns a copy of mol after a constrained MMFF94 relaxation (UFF if
  // MMFF94 setup fails): atoms with index < n_original_atoms are fixed,
  // the appended H atoms move. Requires the bond connectivity
  // (QMAtom::getBondedPartnerIds()), from which OpenBabel's
  // PerceiveBondOrders() derives bond orders.
  //
  // n_steps is the maximum number of conjugate-gradient steps; the
  // relaxation stops earlier once the energy has converged.
  static QMMolecule RelaxNewAtoms(const QMMolecule& mol, Index n_original_atoms,
                                  Index n_steps = 500);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_FRAGMENTSATURATOR_H
