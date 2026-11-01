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
#include <algorithm>
#include <cmath>
#include <iostream>
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

  // Original atoms keep their indices; new H atoms are appended after
  // them, so atom indices computed before saturation stay valid.
  for (const QMAtom& atom : mol) {
    result.push_back(atom);
  }

  // -1 for original atoms (see SaturationResult).
  std::vector<Index> new_atom_parent_ids(mol.size(), -1);

  double bond_length_bohr = bond_length_angstrom * tools::conv::ang2bohr;
  Index new_index = mol.size();
  for (const QMAtom& atom : mol) {
    if (!atom.hasExternalBond()) {
      continue;
    }
    // getExternalBondDirection() is a unit vector: it is normalized in
    // setExternalBondDirection() and only rotated by SegmentMapper.
    Eigen::Vector3d new_pos =
        atom.getPos() + bond_length_bohr * atom.getExternalBondDirection();
    QMAtom new_h(new_index, "H", new_pos);
    // Record the new H-parent bond on both atoms. RelaxNewAtoms builds
    // the OpenBabel connectivity from getBondedPartnerIds(); without
    // this bond the H would be relaxed as an isolated atom. The parent
    // is updated in result, where it has the same index as in mol.
    new_h.AddBondedPartner(atom.getId());
    result[atom.getId()].AddBondedPartner(new_index);
    result.push_back(new_h);
    new_atom_parent_ids.push_back(atom.getId());
    new_index++;
  }

  return SaturationResult{result, new_atom_parent_ids};
}

QMMolecule FragmentSaturator::RelaxNewAtoms(const QMMolecule& mol,
                                            Index n_original_atoms,
                                            Index n_steps) {
  // Build the OBMol from the known connectivity
  // (QMAtom::getBondedPartnerIds(), set by Md2QmEngine/SegmentMapper).
  // OBMol::PerceiveBondOrders() requires the connectivity as input, and
  // geometry-only bond perception is unreliable for newly added atoms.
  OpenBabel::OBMol obmol;
  obmol.BeginModify();
  tools::Elements elements;
  for (const QMAtom& atom : mol) {
    OpenBabel::OBAtom* obatom = obmol.NewAtom();
    obatom->SetAtomicNum(int(elements.getNucCrg(atom.getElement())));
    // OpenBabel uses Angstrom, xtp uses Bohr.
    Eigen::Vector3d pos_angstrom = atom.getPos() / tools::conv::ang2bohr;
    obatom->SetVector(pos_angstrom.x(), pos_angstrom.y(), pos_angstrom.z());
  }

  // OBMol atom indices are 1-based. Bond order 1 is a placeholder;
  // PerceiveBondOrders() derives the actual orders from connectivity and
  // geometry.
  for (const QMAtom& atom : mol) {
    const Index* partners = atom.getBondedPartnerIds();
    for (Index i = 0; i < QMAtom::kMaxBondedPartners; i++) {
      Index partner_id = partners[i];
      // Bonds are stored on both atoms; add each once, from the lower id.
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
    // MMFF94 atom typing can fail, e.g. on aromatic-ring kekulization
    // (openbabel/openbabel#2567); fall back to the more general UFF.
    pFF = OpenBabel::OBForceField::FindForceField("UFF");
    if (pFF == nullptr || !pFF->Setup(obmol)) {
      throw std::runtime_error(
          "FragmentSaturator::RelaxNewAtoms: could not set up either "
          "MMFF94 or UFF for this fragment.");
    }
  }

  // Fix all original atoms; only the appended H atoms
  // (index >= n_original_atoms) move. Uses OBFFConstraints with
  // Setup(mol, constraints), as documented in OpenBabel's forcefield.cpp.
  OpenBabel::OBFFConstraints constraints;
  for (Index i = 0; i < n_original_atoms; i++) {
    constraints.AddAtomConstraint(int(i) + 1);
  }
  if (!pFF->Setup(obmol, constraints)) {
    throw std::runtime_error(
        "FragmentSaturator::RelaxNewAtoms: could not set up "
        "constraints.");
  }

  // OpenBabel ignores the econv argument of ConjugateGradients() and
  // ConjugateGradientsInitialize() (openbabel/openbabel#1366, #2804), so
  // all n_steps would always run. The steps after convergence only
  // accumulate platform-dependent rounding, which caused differences in
  // the final couplings between platforms. The energy criterion of
  // forcefield.cpp (IsNear(e_n2, e_n1, econv)) is therefore applied here,
  // with OpenBabel's default econv = 1e-6.
  double econv = 1e-6;
  pFF->ConjugateGradientsInitialize(int(n_steps), econv);
  double e_prev = pFF->Energy();
  // Check every 10 steps to limit the cost of Energy() evaluations.
  const Index check_interval = 10;
  bool still_running = true;
  Index step = 0;
  bool converged = false;
  for (; step < n_steps && still_running; step += check_interval) {
    Index steps_this_round = std::min(check_interval, n_steps - step);
    still_running = pFF->ConjugateGradientsTakeNSteps(int(steps_this_round));
    double e_now = pFF->Energy();
    if (std::abs(e_now - e_prev) < econv) {
      converged = true;
      step += steps_this_round;
      break;
    }
    e_prev = e_now;
  }
  // Temporary diagnostic. RelaxNewAtoms is static and has no Logger, so
  // std::cerr is used.
  std::cerr << "[RelaxNewAtoms] took " << step << "/" << n_steps
            << " conjugate-gradient steps ("
            << (converged ? "converged early"
                          : (still_running ? "exhausted full budget"
                                           : "OpenBabel itself stopped"))
            << ")" << std::endl;
  pFF->GetCoordinates(obmol);

  // Same elements and ids; only the free H atoms move (fixed atoms are
  // unchanged up to rounding, not re-verified here).
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
