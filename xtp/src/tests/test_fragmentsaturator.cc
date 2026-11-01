/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE fragmentsaturator_test

// Standard includes
#include <chrono>
#include <cstdio>
#include <cstdlib>

// Third party includes
#include <boost/test/unit_test.hpp>
#include <openbabel/atom.h>
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/fragmentsaturator.h"

using namespace votca::xtp;
using namespace votca;
BOOST_AUTO_TEST_SUITE(fragmentsaturator_test)

BOOST_AUTO_TEST_CASE(saturate_single_external_bond_test) {
  // Hand-checkable geometry; FragmentSaturator works on an already
  // mapped QMMolecule, so no mapping files are needed.
  QMMolecule mol("Test", 0);
  QMAtom atom0(0, "C", Eigen::Vector3d::Zero());
  atom0.setExternalBondDirection(Eigen::Vector3d::UnitX());
  QMAtom atom1(1, "H", Eigen::Vector3d::UnitX());
  QMAtom atom2(2, "H", -Eigen::Vector3d::UnitX());
  mol.push_back(atom0);
  mol.push_back(atom1);
  mol.push_back(atom2);

  double bond_length_angstrom = 1.09;
  FragmentSaturator::SaturationResult saturation_result =
      FragmentSaturator::SaturateExternalBonds(mol, bond_length_angstrom);
  const QMMolecule& result = saturation_result.mol;

  // The input is not modified.
  BOOST_CHECK_EQUAL(mol.size(), 3);

  // Exactly one new atom appended (only atom0 has an external bond).
  BOOST_CHECK_EQUAL(result.size(), 4);

  // Original atoms are copied unchanged, at their original indices.
  BOOST_CHECK_EQUAL(result[0].getElement(), "C");
  BOOST_CHECK_EQUAL(result[0].getId(), 0);
  BOOST_CHECK_EQUAL(result[0].getPos().isApprox(Eigen::Vector3d::Zero()), true);

  BOOST_CHECK_EQUAL(result[1].getElement(), "H");
  BOOST_CHECK_EQUAL(result[1].getId(), 1);
  BOOST_CHECK_EQUAL(result[1].getPos().isApprox(Eigen::Vector3d::UnitX()),
                    true);

  BOOST_CHECK_EQUAL(result[2].getElement(), "H");
  BOOST_CHECK_EQUAL(result[2].getId(), 2);
  BOOST_CHECK_EQUAL(result[2].getPos().isApprox(-Eigen::Vector3d::UnitX()),
                    true);

  // New H appended at index 3, at the bond length (in Bohr) along the
  // external-bond direction +X of atom0.
  BOOST_CHECK_EQUAL(result[3].getElement(), "H");
  BOOST_CHECK_EQUAL(result[3].getId(), 3);
  Eigen::Vector3d expected_new_pos =
      Eigen::Vector3d::Zero() +
      (bond_length_angstrom * tools::conv::ang2bohr) * Eigen::Vector3d::UnitX();
  BOOST_CHECK_EQUAL(result[3].getPos().isApprox(expected_new_pos, 1e-8), true);
  if (!result[3].getPos().isApprox(expected_new_pos, 1e-8)) {
    std::cout << "new H position" << std::endl;
    std::cout << result[3].getPos() << std::endl;
    std::cout << "expected" << std::endl;
    std::cout << expected_new_pos << std::endl;
  }

  // The new bond is recorded on both atoms, as required for the
  // OpenBabel connectivity in RelaxNewAtoms.
  BOOST_CHECK_EQUAL(result[3].getBondedPartnerIds()[0], 0);
  BOOST_CHECK_EQUAL(result[0].getBondedPartnerIds()[0], 3);

  // new_atom_parent_ids: -1 for original atoms, parent id 0 for the H.
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids.size(), 4);
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids[0], -1);
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids[1], -1);
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids[2], -1);
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids[3], 0);
}

BOOST_AUTO_TEST_CASE(saturate_no_external_bonds_test) {
  // Without external bonds no atoms are added.
  QMMolecule mol("Test", 0);
  QMAtom atom0(0, "C", Eigen::Vector3d::Zero());
  QMAtom atom1(1, "H", Eigen::Vector3d::UnitX());
  mol.push_back(atom0);
  mol.push_back(atom1);

  FragmentSaturator::SaturationResult saturation_result =
      FragmentSaturator::SaturateExternalBonds(mol);
  BOOST_CHECK_EQUAL(saturation_result.mol.size(), 2);
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids.size(), 2);
}

BOOST_AUTO_TEST_CASE(relax_new_atom_test) {
  // The MMFF94 result is not checked against exact values, only for
  // physically reasonable behavior: the fixed atoms stay in place and a
  // bad H position is corrected.
  //
  // Methane geometry as in DataFiles/segmentmapper/molecule.xyz, with
  // the fourth H at half the C-H bond length (0.545 instead of about
  // 1.09 Angstrom).
  double ang2bohr = tools::conv::ang2bohr;
  Eigen::Vector3d c_pos = Eigen::Vector3d::Zero();
  Eigen::Vector3d h1_pos = Eigen::Vector3d(0.629, 0.629, 0.629) * ang2bohr;
  Eigen::Vector3d h2_pos = Eigen::Vector3d(-0.629, -0.629, 0.629) * ang2bohr;
  Eigen::Vector3d h3_pos = Eigen::Vector3d(0.629, -0.629, -0.629) * ang2bohr;
  Eigen::Vector3d h4_real_direction =
      Eigen::Vector3d(-0.629, 0.629, -0.629).normalized();
  // Deliberately bad: half the real ~1.09 Angstrom C-H bond length.
  Eigen::Vector3d h4_bad_pos = h4_real_direction * 0.545 * ang2bohr;

  QMMolecule mol("Test", 0);
  QMAtom c_atom(0, "C", c_pos);
  QMAtom h1(1, "H", h1_pos);
  QMAtom h2(2, "H", h2_pos);
  QMAtom h3(3, "H", h3_pos);
  QMAtom h4(4, "H", h4_bad_pos);

  // C bonded to all four H, stored on both atoms as in
  // Md2QmEngine::map().
  c_atom.AddBondedPartner(1);
  c_atom.AddBondedPartner(2);
  c_atom.AddBondedPartner(3);
  c_atom.AddBondedPartner(4);
  h1.AddBondedPartner(0);
  h2.AddBondedPartner(0);
  h3.AddBondedPartner(0);
  h4.AddBondedPartner(0);

  mol.push_back(c_atom);
  mol.push_back(h1);
  mol.push_back(h2);
  mol.push_back(h3);
  mol.push_back(h4);

  // Atom 4 (h4) is the only "new" atom -- indices 0-3 are fixed.
  QMMolecule result = FragmentSaturator::RelaxNewAtoms(mol, 4);

  BOOST_CHECK_EQUAL(result.size(), 5);

  // Fixed atoms may move slightly in the minimizer; check an absolute
  // distance of 0.01 Angstrom. Eigen's isApprox() is relative and would
  // require exact equality for c_pos = 0.
  double loose_tol_bohr = 0.01 * ang2bohr;
  BOOST_CHECK_LT((result[0].getPos() - c_pos).norm(), loose_tol_bohr);
  BOOST_CHECK_LT((result[1].getPos() - h1_pos).norm(), loose_tol_bohr);
  BOOST_CHECK_LT((result[2].getPos() - h2_pos).norm(), loose_tol_bohr);
  BOOST_CHECK_LT((result[3].getPos() - h3_pos).norm(), loose_tol_bohr);

  // The new H moved away from its starting position.
  double moved_distance_bohr = (result[4].getPos() - h4_bad_pos).norm();
  BOOST_CHECK_GT(moved_distance_bohr, 0.05 * ang2bohr);

  // Relaxed C-H bond length in a chemically reasonable range.
  double relaxed_bond_length_angstrom =
      (result[4].getPos() - result[0].getPos()).norm() / ang2bohr;
  BOOST_CHECK_GT(relaxed_bond_length_angstrom, 0.9);
  BOOST_CHECK_LT(relaxed_bond_length_angstrom, 1.3);
  if (relaxed_bond_length_angstrom <= 0.9 ||
      relaxed_bond_length_angstrom >= 1.3) {
    std::cout << "relaxed C-H bond length (Angstrom): "
              << relaxed_bond_length_angstrom << std::endl;
  }
}

// Diagnostic timing test, disabled by default; run with
// --run_test=fragmentsaturator_test/relax_new_atom_timing_breakdown_test.
// Times the individual steps of RelaxNewAtoms for the setup of
// relax_new_atom_test (which once took about 42 s on some machines).
BOOST_AUTO_TEST_CASE(relax_new_atom_timing_breakdown_test,
                     *boost::unit_test::disabled()) {
  double ang2bohr = tools::conv::ang2bohr;
  Eigen::Vector3d c_pos = Eigen::Vector3d::Zero();
  Eigen::Vector3d h1_pos = Eigen::Vector3d(0.629, 0.629, 0.629) * ang2bohr;
  Eigen::Vector3d h2_pos = Eigen::Vector3d(-0.629, -0.629, 0.629) * ang2bohr;
  Eigen::Vector3d h3_pos = Eigen::Vector3d(0.629, -0.629, -0.629) * ang2bohr;
  Eigen::Vector3d h4_real_direction =
      Eigen::Vector3d(-0.629, 0.629, -0.629).normalized();
  Eigen::Vector3d h4_bad_pos = h4_real_direction * 0.545 * ang2bohr;

  QMMolecule mol("Test", 0);
  QMAtom c_atom(0, "C", c_pos);
  QMAtom h1(1, "H", h1_pos);
  QMAtom h2(2, "H", h2_pos);
  QMAtom h3(3, "H", h3_pos);
  QMAtom h4(4, "H", h4_bad_pos);
  c_atom.AddBondedPartner(1);
  c_atom.AddBondedPartner(2);
  c_atom.AddBondedPartner(3);
  c_atom.AddBondedPartner(4);
  h1.AddBondedPartner(0);
  h2.AddBondedPartner(0);
  h3.AddBondedPartner(0);
  h4.AddBondedPartner(0);
  mol.push_back(c_atom);
  mol.push_back(h1);
  mol.push_back(h2);
  mol.push_back(h3);
  mol.push_back(h4);

  using clock = std::chrono::steady_clock;
  auto ms = [](clock::duration d) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(d).count();
  };

  auto t_total_start = clock::now();

  // Step 1: build the OBMol (atoms + bonds) -- mirrors
  // RelaxNewAtoms's own first loop exactly.
  auto t0 = clock::now();
  OpenBabel::OBMol obmol;
  obmol.BeginModify();
  tools::Elements elements;
  for (const QMAtom& atom : mol) {
    OpenBabel::OBAtom* obatom = obmol.NewAtom();
    obatom->SetAtomicNum(int(elements.getNucCrg(atom.getElement())));
    Eigen::Vector3d pos_angstrom = atom.getPos() / ang2bohr;
    obatom->SetVector(pos_angstrom.x(), pos_angstrom.y(), pos_angstrom.z());
  }
  for (const QMAtom& atom : mol) {
    const Index* partners = atom.getBondedPartnerIds();
    for (Index i = 0; i < QMAtom::kMaxBondedPartners; i++) {
      Index partner_id = partners[i];
      if (partner_id == -1 || partner_id <= atom.getId()) {
        continue;
      }
      obmol.AddBond(int(atom.getId()) + 1, int(partner_id) + 1, 1);
    }
  }
  obmol.EndModify();
  auto t1 = clock::now();
  std::cout << "[timing] OBMol construction: " << ms(t1 - t0) << " ms"
            << std::endl;

  // Step 2: PerceiveBondOrders.
  obmol.PerceiveBondOrders();
  auto t2 = clock::now();
  std::cout << "[timing] PerceiveBondOrders: " << ms(t2 - t1) << " ms"
            << std::endl;

  // Step 3: FindForceField + first Setup (no constraints).
  OpenBabel::OBForceField* pFF =
      OpenBabel::OBForceField::FindForceField("MMFF94");
  bool first_setup_ok = (pFF != nullptr) && pFF->Setup(obmol);
  auto t3 = clock::now();
  std::cout << "[timing] FindForceField+Setup(no constraints): " << ms(t3 - t2)
            << " ms (ok=" << first_setup_ok << ")" << std::endl;
  BOOST_REQUIRE(first_setup_ok);

  // Step 4: constrained Setup.
  OpenBabel::OBFFConstraints constraints;
  for (Index i = 0; i < 4; i++) {
    constraints.AddAtomConstraint(int(i) + 1);
  }
  bool constrained_setup_ok = pFF->Setup(obmol, constraints);
  auto t4 = clock::now();
  std::cout << "[timing] Constrained Setup: " << ms(t4 - t3)
            << " ms (ok=" << constrained_setup_ok << ")" << std::endl;
  BOOST_REQUIRE(constrained_setup_ok);

  // Step 5: ConjugateGradients with the default n_steps = 500.
  // RelaxNewAtoms stops early on convergence instead, so this is the
  // worst-case timing.
  pFF->ConjugateGradients(500);
  auto t5 = clock::now();
  std::cout << "[timing] ConjugateGradients(500): " << ms(t5 - t4) << " ms"
            << std::endl;

  // Step 6: GetCoordinates.
  pFF->GetCoordinates(obmol);
  auto t6 = clock::now();
  std::cout << "[timing] GetCoordinates: " << ms(t6 - t5) << " ms" << std::endl;

  std::cout << "[timing] TOTAL: " << ms(t6 - t_total_start) << " ms"
            << std::endl;
}
BOOST_AUTO_TEST_SUITE_END()
