/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE test_hdf5
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <cassert>
#include <votca/xtp/checkpoint.h>
#include <votca/xtp/checkpointreader.h>
#include <votca/xtp/checkpointwriter.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmatom.h>

BOOST_AUTO_TEST_SUITE(test_hdf5)
using namespace votca::xtp;
BOOST_AUTO_TEST_CASE(checkpoint_file_test) {

  int basisSetSize = 17;
  int occupiedLevels = 4;
  int unoccupiedLevels = 13;
  int numElectrons = 12;

  Eigen::VectorXd moeTest = Eigen::VectorXd::Random(17);
  Eigen::MatrixXd mocTest = Eigen::MatrixXd::Random(17, 17);

  QMMolecule atoms = QMMolecule(" ", 0);
  for (int i = 0; i < 10; ++i) {
    atoms.push_back(QMAtom(0, "O", Eigen::Vector3d::Random()));
    atoms.push_back(QMAtom(25, "O", Eigen::Vector3d::Random()));
    atoms.push_back(QMAtom(32, "O", Eigen::Vector3d::Random()));
    atoms.push_back(QMAtom(100, "O", Eigen::Vector3d::Random()));
    atoms.push_back(QMAtom(2, "Si", Eigen::Vector3d::Random()));
    atoms.push_back(QMAtom(3145, "N", Eigen::Vector3d::Random()));
  }

  double qmEnergy = -2.1025e-3;

  std::string qmPackage = "NOPE";
  double selfEnergy = 3.14159e23;

  std::string dftBasis = "AWESOME basis*,, 2/.8";
  std::string auxBasis = "cos(theta) = pretty okay basis";

  int rpaMin = '?';
  int rpaMax = 1e3;

  int bseVmin = -6019386;
  int bseCmax = 42;

  double scaHfx = 3.14159;

  bool useTDA = true;

  Eigen::MatrixXd vxcTest = Eigen::MatrixXd::Random(200, 200);
  std::string someECP = "aye aye Cap'n";

  Eigen::MatrixXd QPpertEnergiesTest = Eigen::MatrixXd::Random(31, 42);
  Eigen::MatrixXd QPdiagEnergiesTest = Eigen::VectorXd::Random(21);
  Eigen::MatrixXd QPdiagCoefficientsTest = Eigen::MatrixXd::Identity(42, 42);

  Eigen::VectorXd BSESingletEnergiesTest = Eigen::VectorXd::Random(25);
  Eigen::MatrixXd BSESingletCoefficientsTest = Eigen::MatrixXd::Random(25, 38);
  Eigen::MatrixXd BSESingletCoefficientsARTest =
      Eigen::MatrixXd::Random(42, 42);

  Eigen::VectorXd BSETripletEnergiesTest = Eigen::VectorXd::Random(33);
  Eigen::MatrixXd BSETripletCoefficientsTest = Eigen::MatrixXd::Random(33, 31);

  {
    // Write orbitals
    Orbitals orbWrite;

    orbWrite.setBasisSetSize(basisSetSize);
    orbWrite.setNumberOfOccupiedLevels(occupiedLevels);
    orbWrite.setNumberOfAlphaElectrons(numElectrons);
    orbWrite.MOs().eigenvalues() = moeTest;
    orbWrite.MOs().eigenvectors() = mocTest;

    orbWrite.QMAtoms() = atoms;
    orbWrite.setQMEnergy(qmEnergy);
    orbWrite.setQMpackage(qmPackage);
    orbWrite.setSelfEnergy(selfEnergy);
    orbWrite.setDFTbasisName(dftBasis);
    orbWrite.setAuxbasisName(auxBasis);
    orbWrite.setRPAindices(rpaMin, rpaMax);
    // no need to write qpmin, qpmax
    orbWrite.setBSEindices(bseVmin, bseCmax);
    orbWrite.setScaHFX(scaHfx);
    orbWrite.setTDAApprox(useTDA);
    orbWrite.setECPName(someECP);
    orbWrite.QPpertEnergies() = QPpertEnergiesTest;
    orbWrite.QPdiag().eigenvalues() = QPdiagEnergiesTest;
    orbWrite.QPdiag().eigenvectors() = QPdiagCoefficientsTest;
    orbWrite.BSESinglets().eigenvalues() = BSESingletEnergiesTest;
    orbWrite.BSESinglets().eigenvectors() = BSESingletCoefficientsTest;
    orbWrite.BSESinglets().eigenvectors2() = BSESingletCoefficientsARTest;
    orbWrite.BSETriplets().eigenvalues() = BSETripletEnergiesTest;
    orbWrite.BSETriplets().eigenvectors() = BSETripletCoefficientsTest;

    orbWrite.WriteToCpt("xtp_testing.hdf5");
  }
  // Read Orbitals
  Orbitals orbRead;
  orbRead.ReadFromCpt("xtp_testing.hdf5");

  double tol = 1e-6;

  // Test the read values
  BOOST_CHECK_EQUAL(orbRead.getBasisSetSize(), basisSetSize);
  BOOST_CHECK_EQUAL(orbRead.getBasisSetSize(),
                    occupiedLevels + unoccupiedLevels);
  BOOST_CHECK_EQUAL(orbRead.getNumberOfAlphaElectrons(), numElectrons);
  BOOST_CHECK(orbRead.MOs().eigenvalues().isApprox(moeTest, tol));

  BOOST_CHECK(orbRead.MOs().eigenvectors().isApprox(mocTest, tol));
  BOOST_CHECK_CLOSE(orbRead.getDFTTotalEnergy(), qmEnergy, tol);
  BOOST_CHECK_EQUAL(orbRead.getQMpackage(), qmPackage);
  BOOST_CHECK_CLOSE(orbRead.getSelfEnergy(), selfEnergy, tol);
  BOOST_CHECK_EQUAL(orbRead.getDFTbasisName(), dftBasis);
  BOOST_CHECK_EQUAL(orbRead.getAuxbasisName(), auxBasis);
  BOOST_CHECK_EQUAL(orbRead.getRPAmin(), rpaMin);
  BOOST_CHECK_EQUAL(orbRead.getRPAmax(), rpaMax);

  BOOST_CHECK_EQUAL(orbRead.getBSEvmin(), bseVmin);
  BOOST_CHECK_EQUAL(orbRead.getBSEcmax(), bseCmax);

  BOOST_CHECK_CLOSE(orbRead.getScaHFX(), scaHfx, tol);
  BOOST_CHECK_EQUAL(orbRead.getTDAApprox(), useTDA);
  BOOST_CHECK_EQUAL(orbRead.getECPName(), someECP);
  BOOST_CHECK(orbRead.QPpertEnergies().isApprox(QPpertEnergiesTest, tol));
  BOOST_CHECK(orbRead.QPdiag().eigenvalues().isApprox(QPdiagEnergiesTest, tol));
  BOOST_CHECK(orbRead.QPdiag().eigenvectors().isApprox(QPdiagCoefficientsTest));
  BOOST_CHECK(orbRead.BSESinglets().eigenvalues().isApprox(
      BSESingletEnergiesTest, tol));
  BOOST_CHECK(orbRead.BSESinglets().eigenvectors().isApprox(
      BSESingletCoefficientsTest, tol));
  BOOST_CHECK(orbRead.BSESinglets().eigenvectors2().isApprox(
      BSESingletCoefficientsARTest, tol));
  BOOST_CHECK(orbRead.BSETriplets().eigenvalues().isApprox(
      BSETripletEnergiesTest, tol));
  BOOST_CHECK(orbRead.BSETriplets().eigenvectors().isApprox(
      BSETripletCoefficientsTest, tol));

  BOOST_REQUIRE_EQUAL(orbRead.QMAtoms().size(), atoms.size());

  for (int i = 0; i < atoms.size(); ++i) {
    const auto& atomRead = orbRead.QMAtoms()[i];
    const auto& atomTest = atoms[i];
    BOOST_CHECK_EQUAL(atomRead.getId(), atomTest.getId());
    BOOST_CHECK(atomRead.getPos().isApprox(atomTest.getPos(), tol));
    BOOST_CHECK_EQUAL(atomRead.getNuccharge(), atomTest.getNuccharge());
    BOOST_CHECK_EQUAL(atomRead.getElement(), atomTest.getElement());
  }
}

BOOST_AUTO_TEST_CASE(open_file_error) {
  BOOST_REQUIRE_THROW(
      CheckpointFile cpf("/bin/mr/root/man.pls", CheckpointAccessLevel::READ),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(checkpoint_open_non_existing_loc) {
  CheckpointFile cpf("testin_yo.ab", CheckpointAccessLevel::MODIFY);
  BOOST_REQUIRE_THROW(CheckpointReader r = cpf.getReader("/some/bulshit"),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(read_non_exisiting_matrix) {

  CheckpointFile cpf("xtp_testing.hdf5", CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader("/QMdata");

  Eigen::MatrixXd someMatrix;

  BOOST_REQUIRE_THROW(r(someMatrix, "someMatrix012'5915.jb"),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(read_non_existing_scalar) {
  CheckpointFile cpf("xtp_testing.hdf5", CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader("/QMdata");

  float someThing = 0;
  BOOST_REQUIRE_THROW(r(someThing, "someThing"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(read_vector_strings) {
  CheckpointFile cpf("xtp_vector_string.hdf5");
  CheckpointWriter w = cpf.getWriter();

  std::vector<std::string> test_vec = {
      "a", "b", "advhsbsavc", "dumuunnadnjdsahjads", "DASASDFAFNADDH blasndd"};
  w(test_vec, "vec");

  std::vector<std::string> test_vec2;
  CheckpointReader r = cpf.getReader();
  r(test_vec2, "vec");

  BOOST_CHECK_EQUAL(test_vec.size(), test_vec2.size());

  for (unsigned i = 0; i > test_vec.size(); i++) {
    BOOST_CHECK_EQUAL(test_vec[i], test_vec2[i]);
  }
}

BOOST_AUTO_TEST_CASE(staticsegment) {

  CheckpointFile cpf("xtp_staticsegment.hdf5");
  CheckpointWriter w = cpf.getWriter();
  StaticSegment seg = StaticSegment("test", 0);
  for (int i = 0; i < 10; ++i) {
    seg.push_back(StaticSite(0, "O", Eigen::Vector3d::Random()));
    seg.push_back(StaticSite(25, "O", Eigen::Vector3d::Random()));
    seg.push_back(StaticSite(32, "O", Eigen::Vector3d::Random()));
    seg.push_back(StaticSite(100, "O", Eigen::Vector3d::Random()));
    seg.push_back(StaticSite(2, "Si", Eigen::Vector3d::Random()));
    seg.push_back(StaticSite(3145, "N", Eigen::Vector3d::Random()));
  }

  seg.WriteToCpt(w);
  CheckpointReader r = cpf.getReader();
  StaticSegment seg2 = StaticSegment(r);

  BOOST_REQUIRE_EQUAL(seg2.size(), seg.size());
  for (int i = 0; i < seg.size(); ++i) {
    const auto& atomRead = seg2[i];
    const auto& atomTest = seg[i];
    BOOST_CHECK_EQUAL(atomRead.getId(), atomTest.getId());
    BOOST_CHECK(atomRead.getPos().isApprox(atomTest.getPos(), 1e-7));
    BOOST_CHECK_EQUAL(atomRead.getCharge(), atomTest.getCharge());
    BOOST_CHECK_EQUAL(atomRead.getElement(), atomTest.getElement());
  }
}

BOOST_AUTO_TEST_SUITE_END()
