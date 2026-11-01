/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE qmmolecule_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/qmmolecule.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(qmmolecule_test)

BOOST_AUTO_TEST_CASE(constructors_test) { QMMolecule seg("seg1", 1); }

BOOST_AUTO_TEST_CASE(load_xyz_test) {
  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  QMMolecule seg("seg1", 1);
  seg.LoadFromFile("molecule.xyz");

  auto extension = seg.CalcSpatialMinMax();
  Eigen::Vector3d max(0.629118, 0.629118, 0.629118);
  max *= votca::tools::conv::ang2bohr;
  Eigen::Vector3d min = -max;

  bool check_min = extension.first.isApprox(min, 0.00001);

  if (!check_min) {
    std::cout << "ref" << std::endl;
    std::cout << min << std::endl;
    std::cout << "result" << std::endl;
    std::cout << extension.first << std::endl;
  }
  BOOST_CHECK_EQUAL(check_min, true);
  bool check_max = extension.second.isApprox(max, 0.00001);

  if (!check_max) {
    std::cout << "ref" << std::endl;
    std::cout << max << std::endl;
    std::cout << "result" << std::endl;
    std::cout << extension.second << std::endl;
  }
  BOOST_CHECK_EQUAL(check_max, true);

  seg.WriteXYZ("moltest.xyz", "heloo");
}

BOOST_AUTO_TEST_CASE(readwritehdf) {
  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  QMMolecule seg("seg1", 1);
  seg.LoadFromFile("molecule.xyz");

  CheckpointFile ff("qmmolecule.hdf5");
  CheckpointWriter ww = ff.getWriter();
  seg.WriteToCpt(ww);

  CheckpointReader rr = ff.getReader();
  QMMolecule seg2(rr);
}

BOOST_AUTO_TEST_CASE(
    addcontainer_preserves_external_bond_and_offsets_ids_test) {
  // AddContainer must keep the external-bond data and bonded partners of
  // the merged atoms (they were lost when atoms were reconstructed from
  // index, element and position).
  QMMolecule mol_a("A", 0);
  QMAtom a0(0, "C", Eigen::Vector3d::Zero());
  mol_a.push_back(a0);

  QMMolecule mol_b("B", 1);
  QMAtom b0(0, "S", Eigen::Vector3d::UnitX());
  b0.setExternalBondDirection(Eigen::Vector3d::UnitY());
  b0.setExternalBondPartnerSegmentId(42);
  b0.AddBondedPartner(1);
  QMAtom b1(1, "H", 2 * Eigen::Vector3d::UnitX());
  b1.AddBondedPartner(0);
  mol_b.push_back(b0);
  mol_b.push_back(b1);

  mol_a.AddContainer(mol_b);

  // IDs of the merged atoms are offset by the original size of mol_a (1).
  BOOST_CHECK_EQUAL(mol_a.size(), 3);
  BOOST_CHECK_EQUAL(mol_a[1].getId(), 1);
  BOOST_CHECK_EQUAL(mol_a[2].getId(), 2);

  // External-bond data does not refer to atoms of the container and is
  // kept unchanged.
  BOOST_CHECK_EQUAL(mol_a[1].hasExternalBond(), true);
  BOOST_CHECK_EQUAL(
      mol_a[1].getExternalBondDirection().isApprox(Eigen::Vector3d::UnitY()),
      true);
  BOOST_CHECK_EQUAL(mol_a[1].getExternalBondPartnerSegmentId(), 42);
  BOOST_CHECK_EQUAL(mol_a[2].hasExternalBond(), false);

  // Bonded-partner IDs are local to mol_b and are offset by 1:
  // b0: 1 -> 2, b1: 0 -> 1.
  BOOST_CHECK_EQUAL(mol_a[1].getBondedPartnerIds()[0], 2);
  BOOST_CHECK_EQUAL(mol_a[2].getBondedPartnerIds()[0], 1);
}

BOOST_AUTO_TEST_SUITE_END()
