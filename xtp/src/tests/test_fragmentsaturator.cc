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
  // Simple, direct, hand-verifiable geometry -- deliberately NOT
  // reusing test_segmentmapper.cc's own, more complex, mapping-file-
  // based setup at all, since FragmentSaturator operates directly on
  // an already-mapped QMMolecule and needs none of that machinery.
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

  // mol itself, the original input, must remain entirely unchanged --
  // SaturateExternalBonds must never modify its own input at all.
  BOOST_CHECK_EQUAL(mol.size(), 3);

  // Exactly one new atom appended (only atom0 has an external bond).
  BOOST_CHECK_EQUAL(result.size(), 4);

  // The three original atoms are copied unchanged, at their own,
  // original indices.
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

  // The new H atom itself: appended strictly after all original
  // atoms (index 3), placed along atom0's own, already-normalized
  // external-bond direction (+X), at the given bond length, converted
  // to Bohr the same way SaturateExternalBonds itself does.
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

  // Real, direct, new coverage for a real, direct bug fix: the new H
  // atom's own bonded_partner_ids_ must record its own real parent
  // (atom0), and atom0's own copy within result must, symmetrically,
  // record the new H as one of its own real partners too -- both
  // needed for RelaxNewAtoms's own OpenBabel-based connectivity to
  // actually know these two atoms are bonded at all.
  BOOST_CHECK_EQUAL(result[3].getBondedPartnerIds()[0], 0);
  BOOST_CHECK_EQUAL(result[0].getBondedPartnerIds()[0], 3);

  // new_atom_parent_ids: -1 for every original atom (0-2), and atom0's
  // own id (0) -- the atom the new H is saturating -- for the new H
  // itself (index 3).
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids.size(), 4);
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids[0], -1);
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids[1], -1);
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids[2], -1);
  BOOST_CHECK_EQUAL(saturation_result.new_atom_parent_ids[3], 0);
}

BOOST_AUTO_TEST_CASE(saturate_no_external_bonds_test) {
  // Direct, simple sanity check: if no atom has an external bond at
  // all, the result must be identical to the input (in size -- no new
  // atoms appended at all).
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
  // Sanity-check style test, deliberately not an exact-value one --
  // unlike this file's own earlier tests, the actual, real MMFF94-
  // relaxed position cannot be hand-predicted or independently
  // reproduced here at all (no local OpenBabel available to check
  // against, the way earlier tests in this session could reproduce
  // xtp's own math directly in Python). Checks the real, physically
  // meaningful properties instead: does relaxation actually move a
  // deliberately bad starting position toward something chemically
  // reasonable, while leaving the fixed, original atoms alone.
  //
  // Reuses the same known, tetrahedral methane geometry used
  // elsewhere this session (test_segmentmapper.cc's own
  // DataFiles/segmentmapper/molecule.xyz) -- C at the origin, three
  // real H's at their own, real tetrahedral positions -- but the
  // fourth H is deliberately placed at HALF its own, real distance
  // from the carbon (~0.545 Angstrom instead of the real ~1.09), a
  // clearly, deliberately too-short, chemically wrong starting bond
  // length -- exactly the kind of case relaxation is meant to fix.
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

  // Real, known connectivity: C bonded to all four H's (both
  // directions, matching Md2QmEngine::map()'s own, established,
  // symmetric convention).
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

  // The four originally-fixed atoms come back essentially unchanged
  // -- a real force-field minimizer may move a "fixed" atom by a
  // tiny, real numerical amount depending on its own internal
  // implementation, so this uses a loose, direct, absolute-distance
  // check (0.01 Angstrom) rather than exact equality. Deliberately
  // NOT Eigen's own isApprox() here -- that uses a RELATIVE
  // tolerance (|a-b| <= tol*min(|a|,|b|)), which breaks down against
  // c_pos specifically, exactly Zero() -- a relative tolerance
  // against a zero vector effectively demands exact equality, not
  // the loose check actually intended here.
  double loose_tol_bohr = 0.01 * ang2bohr;
  BOOST_CHECK_LT((result[0].getPos() - c_pos).norm(), loose_tol_bohr);
  BOOST_CHECK_LT((result[1].getPos() - h1_pos).norm(), loose_tol_bohr);
  BOOST_CHECK_LT((result[2].getPos() - h2_pos).norm(), loose_tol_bohr);
  BOOST_CHECK_LT((result[3].getPos() - h3_pos).norm(), loose_tol_bohr);

  // The new H (atom 4) actually moved from its own, deliberately bad
  // starting position.
  double moved_distance_bohr = (result[4].getPos() - h4_bad_pos).norm();
  BOOST_CHECK_GT(moved_distance_bohr, 0.05 * ang2bohr);

  // The resulting C-H bond length for the relaxed atom is now in a
  // chemically reasonable range -- not pinned to a single, exact
  // expected value at all, given no independent way to reproduce
  // OpenBabel's own, real MMFF94 result here.
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

// Direct, temporary diagnostic test -- NOT part of the normal test
// suite at all (disabled by default: run explicitly via
// --run_test=fragmentsaturator_test/relax_new_atom_timing_breakdown_test).
// Added directly in response to a real, observed, surprising ~42s
// runtime for relax_new_atom_test above, on the user's own real
// machine -- re-implements RelaxNewAtoms's own steps here, separately,
// each directly timed, since I have no way to run a real profiler on
// the user's own machine myself. Uses the exact same geometry/setup
// as relax_new_atom_test above, so the timings here are directly
// comparable to that real, already-observed 42s figure.
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

  // Step 5: ConjugateGradients (the same n_steps=500 default
  // RelaxNewAtoms itself uses). Real, direct note: RelaxNewAtoms's
  // own real production code no longer calls this single-shot
  // ConjugateGradients() directly at all (replaced with a manual,
  // early-stopping convergence loop, in direct response to a real,
  // confirmed cross-platform coupling discrepancy this session) -- so
  // this specific step no longer exactly mirrors production. Left
  // as-is deliberately: this is still a genuinely useful worst-case
  // timing baseline (the full 500 steps, with no early stopping at
  // all), which is exactly what this diagnostic test is for.
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
