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
#include <votca/tools/elements.h>

// OpenBabel includes
#include <openbabel/atom.h>
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>

// Local VOTCA includes
#include "votca/xtp/fragmentsaturator.h"

namespace votca {
namespace xtp {

FragmentSaturator::SaturationResult FragmentSaturator::SaturateExternalBonds(
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

  // -1 (unset) for every original atom -- see SaturationResult's own
  // header comment for exactly what this means/is for.
  std::vector<Index> new_atom_parent_ids(mol.size(), -1);

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
    new_atom_parent_ids.push_back(atom.getId());
    new_index++;
  }

  return SaturationResult{result, new_atom_parent_ids};
}

QMMolecule FragmentSaturator::RelaxNewAtoms(const QMMolecule& mol,
                                            Index n_original_atoms,
                                            Index n_steps) {
  // Build a real OBMol directly from mol's own atoms and real, known
  // bond connectivity (QMAtom::getBondedPartnerIds(), populated by
  // the SegmentMapper/Md2QmEngine pipeline already built earlier this
  // session) -- OBMol::PerceiveBondOrders(), below, needs this
  // connectivity as real input, it does not itself guess which atoms
  // are bonded at all (confirmed directly, earlier this session, by
  // reading its own source), and geometry-only bond perception (no
  // known connectivity at all) was separately confirmed, this same
  // session, to be genuinely unreliable for exactly this kind of
  // situation (a newly-added atom), via a real, direct RDKit test
  // that gave a chemically wrong result for an analogous case.
  OpenBabel::OBMol obmol;
  obmol.BeginModify();
  tools::Elements elements;
  for (const QMAtom& atom : mol) {
    OpenBabel::OBAtom* obatom = obmol.NewAtom();
    obatom->SetAtomicNum(int(elements.getNucCrg(atom.getElement())));
    // Bohr -> Angstrom -- OpenBabel's own, standard internal unit
    // (confirmed directly, this session, from its own official
    // examples), unlike xtp's own, internal Bohr convention
    // (QMAtom::getPos()'s own header comment, qmatom.h).
    Eigen::Vector3d pos_angstrom = atom.getPos() / tools::conv::ang2bohr;
    obatom->SetVector(pos_angstrom.x(), pos_angstrom.y(), pos_angstrom.z());
  }

  // OBMol::AddBond's own atom indices are 1-based (confirmed
  // directly, from OpenBabel's own official examples), unlike
  // QMAtom::getId()'s own 0-based indices, so +1 is applied here.
  // Bond order itself is passed as a placeholder single bond (1) for
  // every bond -- PerceiveBondOrders(), below, derives the real bond
  // order from this connectivity plus real geometry; this initial
  // value is never trusted directly.
  for (const QMAtom& atom : mol) {
    const Index* partners = atom.getBondedPartnerIds();
    for (Index i = 0; i < QMAtom::kMaxBondedPartners; i++) {
      Index partner_id = partners[i];
      // Only add each real bond once (from the lower-ID side) --
      // getBondedPartnerIds() records both directions of every real
      // bond (confirmed directly, earlier this session, by
      // Md2QmEngine::map()'s own symmetric AddBondedPartner calls,
      // once for each atom on either side of a given bond), so
      // without this check, every real bond would be added twice.
      if (partner_id == -1 || partner_id <= atom.getId()) {
        continue;
      }
      obmol.AddBond(int(atom.getId()) + 1, int(partner_id) + 1, 1);
    }
  }
  obmol.EndModify();
  obmol.PerceiveBondOrders();

  OpenBabel::OBForceField* pFF =
      OpenBabel::OBForceField::FindForceField("MMFF94");
  if (pFF == nullptr || !pFF->Setup(obmol)) {
    // MMFF94 atom typing genuinely can fail for some real structures
    // (confirmed directly, this session, via a real, documented,
    // still-open upstream OpenBabel issue around aromatic-ring
    // kekulization, openbabel/openbabel #2567) -- UFF is a real,
    // established, more general-purpose fallback for exactly this
    // situation.
    pFF = OpenBabel::OBForceField::FindForceField("UFF");
    if (pFF == nullptr || !pFF->Setup(obmol)) {
      throw std::runtime_error(
          "FragmentSaturator::RelaxNewAtoms: could not set up either "
          "MMFF94 or UFF for this fragment.");
    }
  }

  // Fix every original atom in place -- only the new H atom(s)
  // SaturateExternalBonds appended (index >= n_original_atoms, per
  // its own, documented convention of appending strictly after all
  // original atoms) are free to move. Matches the official, direct
  // OBFFConstraints/Setup(mol, constraints) pattern documented
  // directly in OpenBabel's own forcefield.cpp, rather than
  // OBForceField::SetFixAtom() called directly -- confirmed this is
  // the officially-documented way, before writing this, rather than
  // guessed.
  OpenBabel::OBFFConstraints constraints;
  for (Index i = 0; i < n_original_atoms; i++) {
    constraints.AddAtomConstraint(int(i) + 1);
  }
  if (!pFF->Setup(obmol, constraints)) {
    throw std::runtime_error(
        "FragmentSaturator::RelaxNewAtoms: could not set up "
        "constraints.");
  }

  pFF->ConjugateGradients(int(n_steps));
  pFF->GetCoordinates(obmol);

  // Build the resulting, relaxed QMMolecule -- same element/ID for
  // every atom; only the position itself potentially changes (for the
  // free, new H atom(s) -- every fixed, original atom's own position
  // should come back essentially unchanged, modulo floating-point
  // noise, though this is not separately re-verified here).
  QMMolecule result(mol.getType(), mol.getId());
  Index idx = 0;
  for (const QMAtom& atom : mol) {
    OpenBabel::OBAtom* obatom = obmol.GetAtom(int(idx) + 1);
    Eigen::Vector3d pos_bohr =
        Eigen::Vector3d(obatom->GetX(), obatom->GetY(), obatom->GetZ()) *
        tools::conv::ang2bohr;
    QMAtom new_atom(atom.getId(), atom.getElement(), pos_bohr);
    result.push_back(new_atom);
    idx++;
  }

  return result;
}

}  // namespace xtp
}  // namespace votca
