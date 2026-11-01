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

#define BOOST_TEST_MODULE segmentmapper_test

// Standard includes
#include <cstdio>
#include <cstdlib>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/segmentmapper.h"

using namespace votca::xtp;
using namespace votca;
BOOST_AUTO_TEST_SUITE(segmentmapper_test)

BOOST_AUTO_TEST_CASE(mapping_test) {

  Logger log;
  QMMapper mapper = QMMapper(log);
  mapper.LoadMappingFile(std::string(XTP_TEST_DATA_FOLDER) +
                         "/segmentmapper/ch4.xml");
  Segment seg("Methane", 1);
  Atom atm1(1, "CB", 5, Eigen::Vector3d::Zero(), "C");
  Atom atm2(1, "HB1", 6, Eigen::Vector3d::UnitX(), "H");
  Atom atm3(1, "HB2", 7, Eigen::Vector3d::UnitY(), "H");
  Atom atm4(1, "HB3", 8, -Eigen::Vector3d::UnitX(), "H");
  Atom atm5(1, "HB4", 9, -Eigen::Vector3d::UnitY(), "H");
  seg.push_back(atm1);
  seg.push_back(atm2);
  seg.push_back(atm3);
  seg.push_back(atm4);
  seg.push_back(atm5);

  QMMolecule qmmol = mapper.map(
      seg, std::string(XTP_TEST_DATA_FOLDER) + "/segmentmapper/molecule.xyz");
  std::vector<std::string> name_ref = {"C", "H", "H", "H", "H"};
  std::vector<Index> id_ref = {0, 1, 2, 3, 4};
  std::vector<Eigen::Vector3d> pos_ref;
  Eigen::Vector3d pos1 = {-0.026627, -0.0672429, 8.10559e-19};
  pos_ref.push_back(pos1);
  Eigen::Vector3d pos2 = {2.03254, -0.0672429, 3.24336e-17};
  pos_ref.push_back(pos2);
  Eigen::Vector3d pos3 = {-0.713016, 1.87416, -4.21603e-17};
  pos_ref.push_back(pos3);
  Eigen::Vector3d pos4 = {-1, 0, 0};
  pos_ref.push_back(pos4);
  Eigen::Vector3d pos5 = {0, -1, 0};
  pos_ref.push_back(pos5);
  BOOST_CHECK_EQUAL(qmmol.getId(), 1);
  BOOST_CHECK_EQUAL(qmmol.getType(), "Methane");
  Eigen::Vector3d ref = {0.000140384, 0.000354522, 0};
  BOOST_CHECK_EQUAL(qmmol.getPos().isApprox(ref, 1e-5), true);
  if (!qmmol.getPos().isApprox(ref, 1e-5)) {
    std::cout << "qmmolpos" << std::endl;
    std::cout << qmmol.getPos() << std::endl;
    std::cout << "qmmolref" << std::endl;
    std::cout << ref << std::endl;
  }
  for (Index i = 0; i < qmmol.size(); i++) {
    BOOST_CHECK_EQUAL(qmmol[i].getElement(), name_ref[i]);
    bool pos_equal = qmmol[i].getPos().isApprox(pos_ref[i], 1e-5);
    if (!pos_equal) {
      std::cout << "Atom " << i << std::endl;
      std::cout << "pos" << std::endl;
      std::cout << qmmol[i].getPos() << std::endl;
      std::cout << "ref" << std::endl;
      std::cout << pos_ref[i] << std::endl;
    }
    BOOST_CHECK_EQUAL(qmmol[i].getId(), id_ref[i]);
  }
}

BOOST_AUTO_TEST_CASE(mapping_to_md_test) {

  Logger log;
  QMMapper mapper = QMMapper(log);
  mapper.LoadMappingFile(std::string(XTP_TEST_DATA_FOLDER) +
                         "/segmentmapper/ch4_2.xml");
  Segment seg("Methane", 1);
  Atom atm1(1, "CB", 5, Eigen::Vector3d::Zero(), "C");
  Atom atm2(1, "HB1", 6, Eigen::Vector3d::UnitX(), "H");
  Atom atm3(1, "HB2", 7, Eigen::Vector3d::UnitY(), "H");
  Atom atm4(1, "HB3", 8, -Eigen::Vector3d::UnitX(), "H");
  Atom atm5(1, "HB4", 9, -Eigen::Vector3d::UnitY(), "H");
  seg.push_back(atm1);
  seg.push_back(atm2);
  seg.push_back(atm3);
  seg.push_back(atm4);
  seg.push_back(atm5);

  QMMolecule qmmol = mapper.map(
      seg, std::string(XTP_TEST_DATA_FOLDER) + "/segmentmapper/molecule2.xyz");
  std::vector<std::string> name_ref = {"C", "H", "H", "H", "H"};
  std::vector<Index> id_ref = {0, 1, 2, 3, 4};
  std::vector<Eigen::Vector3d> pos_ref;
  Eigen::Vector3d pos1 = {0, 0, 0};
  pos_ref.push_back(pos1);
  Eigen::Vector3d pos2 = {1, 0, 0};
  pos_ref.push_back(pos2);
  Eigen::Vector3d pos3 = {0, 1, 0};
  pos_ref.push_back(pos3);
  Eigen::Vector3d pos4 = {-1, 0, 0};
  pos_ref.push_back(pos4);
  Eigen::Vector3d pos5 = {0, -1, 0};
  pos_ref.push_back(pos5);
  BOOST_CHECK_EQUAL(qmmol.getId(), 1);
  BOOST_CHECK_EQUAL(qmmol.getType(), "Methane");
  Eigen::Vector3d ref = {0.0, 0.0, 0};
  BOOST_CHECK_EQUAL(qmmol.getPos().isApprox(ref, 1e-2), true);
  if (!qmmol.getPos().isApprox(ref, 1e-5)) {
    std::cout << "qmmolpos" << std::endl;
    std::cout << qmmol.getPos() << std::endl;
    std::cout << "qmmolref" << std::endl;
    std::cout << ref << std::endl;
  }
  for (Index i = 0; i < qmmol.size(); i++) {
    BOOST_CHECK_EQUAL(qmmol[i].getElement(), name_ref[i]);
    bool pos_equal = qmmol[i].getPos().isApprox(pos_ref[i], 1e-5);
    if (!pos_equal) {
      std::cout << "Atom " << i << std::endl;
      std::cout << "pos" << std::endl;
      std::cout << qmmol[i].getPos() << std::endl;
      std::cout << "ref" << std::endl;
      std::cout << pos_ref[i] << std::endl;
    }
    BOOST_CHECK_EQUAL(qmmol[i].getId(), id_ref[i]);
  }
}

BOOST_AUTO_TEST_CASE(mapping_test_no_weights) {

  Logger log;
  QMMapper mapper = QMMapper(log);
  mapper.LoadMappingFile(std::string(XTP_TEST_DATA_FOLDER) +
                         "/segmentmapper/ch4_3.xml");
  Segment seg("Methane", 1);
  Atom atm1(1, "CB", 5, Eigen::Vector3d::Ones(), "C");
  Atom atm2(1, "HB1", 6, Eigen::Vector3d::UnitX() + Eigen::Vector3d::Ones(),
            "H");
  Atom atm3(1, "HB2", 7, Eigen::Vector3d::UnitY() + Eigen::Vector3d::Ones(),
            "H");
  Atom atm4(1, "HB3", 8, -Eigen::Vector3d::UnitX() + Eigen::Vector3d::Ones(),
            "H");
  Atom atm5(1, "HB4", 9, -Eigen::Vector3d::UnitY() + Eigen::Vector3d::Ones(),
            "H");
  seg.push_back(atm1);
  seg.push_back(atm2);
  seg.push_back(atm3);
  seg.push_back(atm4);
  seg.push_back(atm5);

  QMMolecule qmmol = mapper.map(
      seg, std::string(XTP_TEST_DATA_FOLDER) + "/segmentmapper/molecule3.xyz");
  std::vector<std::string> name_ref = {"C", "H", "H", "H", "H"};
  std::vector<Index> id_ref = {0, 1, 2, 3, 4};
  std::vector<Eigen::Vector3d> pos_ref;
  Eigen::Vector3d pos1 = {1.0 - 0.0267876, 1.0 - 0.0676484, 1};
  pos_ref.push_back(pos1);
  Eigen::Vector3d pos2 = {1 + 2.03238, 1.0 - 0.0676484, 1};
  pos_ref.push_back(pos2);
  Eigen::Vector3d pos3 = {1 + -0.713177, 1 + 1.87375, 1};
  pos_ref.push_back(pos3);
  Eigen::Vector3d pos4 = {0, 1, 1};
  pos_ref.push_back(pos4);
  Eigen::Vector3d pos5 = {1, 0, 1};
  pos_ref.push_back(pos5);
  BOOST_CHECK_EQUAL(qmmol.getId(), 1);
  BOOST_CHECK_EQUAL(qmmol.getType(), "Methane");
  Eigen::Vector3d ref = {1, 1, 1};
  BOOST_CHECK_EQUAL(qmmol.getPos().isApprox(ref, 1e-5), true);
  if (!qmmol.getPos().isApprox(ref, 1e-5)) {
    std::cout << "qmmolpos" << std::endl;
    std::cout << qmmol.getPos() << std::endl;
    std::cout << "qmmolref" << std::endl;
    std::cout << ref << std::endl;
  }
  for (Index i = 0; i < qmmol.size(); i++) {
    BOOST_CHECK_EQUAL(qmmol[i].getElement(), name_ref[i]);
    bool pos_equal = qmmol[i].getPos().isApprox(pos_ref[i], 1e-5);
    if (!pos_equal) {
      std::cout << "Atom " << i << std::endl;
      std::cout << "pos" << std::endl;
      std::cout << qmmol[i].getPos() << std::endl;
      std::cout << "ref" << std::endl;
      std::cout << pos_ref[i] << std::endl;
    }
    BOOST_CHECK_EQUAL(qmmol[i].getId(), id_ref[i]);
  }
}

BOOST_AUTO_TEST_CASE(external_bond_direction_transfer_test) {
  // Tests the transfer of the external-bond direction and partner segment
  // ID from Atom to QMAtom in SegmentMapper, using the setup of
  // mapping_test. The template molecule.xyz is a tetrahedral methane while
  // the MD geometry is planar, so the out-of-plane axis is not preserved.
  // The reference direction was obtained by reproducing the rotation of
  // MapMapAtomonMD (rotateMAP2MD = rot_md * rot_map^T) independently in
  // Python.
  Logger log;
  QMMapper mapper = QMMapper(log);
  mapper.LoadMappingFile(std::string(XTP_TEST_DATA_FOLDER) +
                         "/segmentmapper/ch4.xml");
  Segment seg("Methane", 1);
  Atom atm1(1, "CB", 5, Eigen::Vector3d::Zero(), "C");
  Atom atm2(1, "HB1", 6, Eigen::Vector3d::UnitX(), "H");
  Atom atm3(1, "HB2", 7, Eigen::Vector3d::UnitY(), "H");
  Atom atm4(1, "HB3", 8, -Eigen::Vector3d::UnitX(), "H");
  Atom atm5(1, "HB4", 9, -Eigen::Vector3d::UnitY(), "H");

  // Only atm1 ("C") gets an external bond, along +Z. The partner atom ID
  // 999 is a placeholder, distinct from the atom IDs 5-9; its resolution is
  // tested in test_md2qmengine.cc.
  atm1.setExternalBondDirection(Eigen::Vector3d::UnitZ(), 999);
  // Arbitrary partner segment ID, distinct from atom and segment IDs; it
  // must be copied unchanged.
  atm1.setExternalBondPartnerSegmentId(42);

  seg.push_back(atm1);
  seg.push_back(atm2);
  seg.push_back(atm3);
  seg.push_back(atm4);
  seg.push_back(atm5);

  QMMolecule qmmol = mapper.map(
      seg, std::string(XTP_TEST_DATA_FOLDER) + "/segmentmapper/molecule.xyz");

  // Atom 0 of the mapped result is the "C" (see mapping_test).
  BOOST_CHECK_EQUAL(qmmol[0].hasExternalBond(), true);
  Eigen::Vector3d dir = qmmol[0].getExternalBondDirection();

  // Normalization and rotation preserve unit length.
  BOOST_CHECK_CLOSE(dir.norm(), 1.0, 1e-6);

  // Reference value, see comment at the top of this test.
  Eigen::Vector3d dir_ref = {0.57735026918962573, 0.81649658092772615, 0.0};
  BOOST_CHECK_EQUAL(dir.isApprox(dir_ref, 1e-5), true);
  if (!dir.isApprox(dir_ref, 1e-5)) {
    std::cout << "external bond direction" << std::endl;
    std::cout << dir << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << dir_ref << std::endl;
  }

  // The partner segment ID is not geometric and is copied unchanged.
  BOOST_CHECK_EQUAL(qmmol[0].getExternalBondPartnerSegmentId(), 42);

  // The other atoms have no external bond.
  for (Index i = 1; i < qmmol.size(); i++) {
    BOOST_CHECK_EQUAL(qmmol[i].hasExternalBond(), false);
  }
}

BOOST_AUTO_TEST_CASE(bonded_partner_translation_across_fragments_test) {
  // Tests the translation of MD-level bonded-partner IDs to QM-level IDs
  // for partners in different fragments, which requires the second pass in
  // SegmentMapper::map(). ch4.xml splits the segment into three fragments
  // (ALAB_1: atoms 0,1,2; ALAB_2: atom 3; ALAB_3: atom 4), and the C
  // (atom 0) is bonded to atoms in all three. In ch4.xml MD atom k maps to
  // QM atom k (1:CB:0 -> 0:C, 1:HB1:1 -> 1:H, ...), so the expected QM IDs
  // equal the mapping indices.
  Logger log;
  QMMapper mapper = QMMapper(log);
  mapper.LoadMappingFile(std::string(XTP_TEST_DATA_FOLDER) +
                         "/segmentmapper/ch4.xml");
  Segment seg("Methane", 1);
  Atom atm1(1, "CB", 5, Eigen::Vector3d::Zero(), "C");
  Atom atm2(1, "HB1", 6, Eigen::Vector3d::UnitX(), "H");
  Atom atm3(1, "HB2", 7, Eigen::Vector3d::UnitY(), "H");
  Atom atm4(1, "HB3", 8, -Eigen::Vector3d::UnitX(), "H");
  Atom atm5(1, "HB4", 9, -Eigen::Vector3d::UnitY(), "H");

  // Methane connectivity: C bonded to all four H.
  atm1.AddBondedPartner(atm2.getId());
  atm1.AddBondedPartner(atm3.getId());
  atm1.AddBondedPartner(atm4.getId());
  atm1.AddBondedPartner(atm5.getId());
  atm2.AddBondedPartner(atm1.getId());
  atm3.AddBondedPartner(atm1.getId());
  atm4.AddBondedPartner(atm1.getId());
  atm5.AddBondedPartner(atm1.getId());

  seg.push_back(atm1);
  seg.push_back(atm2);
  seg.push_back(atm3);
  seg.push_back(atm4);
  seg.push_back(atm5);

  QMMolecule qmmol = mapper.map(
      seg, std::string(XTP_TEST_DATA_FOLDER) + "/segmentmapper/molecule.xyz");

  // Atom 0 (C): partners 1, 2 (same fragment) and 3, 4 (other fragments),
  // in the order added, since AddBondedPartner fills the first free slot.
  const Index* c_partners = qmmol[0].getBondedPartnerIds();
  BOOST_CHECK_EQUAL(c_partners[0], 1);
  BOOST_CHECK_EQUAL(c_partners[1], 2);
  BOOST_CHECK_EQUAL(c_partners[2], 3);
  BOOST_CHECK_EQUAL(c_partners[3], 4);

  // Each H (atoms 1-4): partner 0, remaining slots unset (-1).
  for (Index i = 1; i < qmmol.size(); i++) {
    const Index* h_partners = qmmol[i].getBondedPartnerIds();
    BOOST_CHECK_EQUAL(h_partners[0], 0);
    BOOST_CHECK_EQUAL(h_partners[1], -1);
    BOOST_CHECK_EQUAL(h_partners[2], -1);
    BOOST_CHECK_EQUAL(h_partners[3], -1);
  }
}
BOOST_AUTO_TEST_SUITE_END()
