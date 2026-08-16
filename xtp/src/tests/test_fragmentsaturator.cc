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
#include <cstdio>
#include <cstdlib>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

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
  QMMolecule result =
      FragmentSaturator::SaturateExternalBonds(mol, bond_length_angstrom);

  // mol itself, the original input, must remain entirely unchanged --
  // SaturateExternalBonds must never modify its own input at all.
  BOOST_CHECK_EQUAL(mol.size(), 3);

  // Exactly one new atom appended (only atom0 has an external bond).
  BOOST_CHECK_EQUAL(result.size(), 4);

  // The three original atoms are copied unchanged, at their own,
  // original indices.
  BOOST_CHECK_EQUAL(result[0].getElement(), "C");
  BOOST_CHECK_EQUAL(result[0].getId(), 0);
  BOOST_CHECK_EQUAL(result[0].getPos().isApprox(Eigen::Vector3d::Zero()),
                    true);

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

  QMMolecule result = FragmentSaturator::SaturateExternalBonds(mol);
  BOOST_CHECK_EQUAL(result.size(), 2);
}
BOOST_AUTO_TEST_SUITE_END()
