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

#define BOOST_TEST_MODULE aomatrix_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/orbitals.h"
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(aomatrix_test)

QMMolecule Methane() {

  QMMolecule mol(" ", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/aomatrix/molecule.xyz");
  return mol;
}

BOOST_AUTO_TEST_CASE(aomatrices_test) {

  QMMolecule mol = Methane();
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  AOOverlap overlap;
  overlap.Fill(aobasis);
  Eigen::MatrixXd overlap_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/overlap_ref.mm");

  bool check_overlap = overlap.Matrix().isApprox(overlap_ref, 0.0001);
  BOOST_CHECK_EQUAL(check_overlap, 1);
  if (!check_overlap) {
    cout << "ref" << endl;
    cout << overlap_ref << endl;
    cout << "result" << endl;
    cout << overlap.Matrix() << endl;
  }

  AOKinetic kinetic;
  kinetic.Fill(aobasis);
  Eigen::MatrixXd kinetic_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/kinetic_ref.mm");

  bool check_kinetic = kinetic.Matrix().isApprox(kinetic_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_kinetic, 1);
  if (!check_kinetic) {
    cout << "ref" << endl;
    cout << kinetic_ref << endl;
    cout << "result" << endl;
    cout << kinetic.Matrix() << endl;
  }

  AOCoulomb coulomb;
  coulomb.Fill(aobasis);
  Eigen::MatrixXd coulomb_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/coulomb_ref.mm");
  bool check_coulomb = coulomb.Matrix().isApprox(coulomb_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_coulomb, 1);
  if (!check_coulomb) {
    cout << "ref" << endl;
    cout << coulomb_ref << endl;
    cout << "result" << endl;
    cout << coulomb.Matrix() << endl;
  }

  Eigen::MatrixXd ps_invSqrt = coulomb.Pseudo_InvSqrt(1e-7);
  Eigen::MatrixXd coulombinvsqrt_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/aomatrix/coulombinvsqrt_ref.mm");

  bool check_coulombinvsqrt = ps_invSqrt.isApprox(coulombinvsqrt_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_coulombinvsqrt, 1);
  if (!check_coulombinvsqrt) {
    cout << "ref" << endl;
    cout << coulombinvsqrt_ref << endl;
    cout << "result" << endl;
    cout << ps_invSqrt << endl;
  }

  Eigen::MatrixXd ps_invSqrtgw = coulomb.Pseudo_InvSqrt_GWBSE(overlap, 1e-7);
  Eigen::MatrixXd coulombinvsqrtgw_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/aomatrix/coulombinvsqrtgw_ref.mm");

  bool check_coulombinvsqrtgw =
      ps_invSqrtgw.isApprox(coulombinvsqrtgw_ref, 0.00001);

  BOOST_CHECK_EQUAL(check_coulombinvsqrtgw, 1);
  if (!check_coulombinvsqrtgw) {
    cout << "ref" << endl;
    cout << coulombinvsqrtgw_ref << endl;
    cout << "result" << endl;
    cout << ps_invSqrtgw << endl;
  }
}

BOOST_AUTO_TEST_CASE(aomatrices_contracted_test) {

  QMMolecule mol("C", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/C.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/contracted.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  AOOverlap overlap;
  overlap.Fill(aobasis);
  Eigen::MatrixXd overlap_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/aomatrix/overlap_ref_contracted.mm");

  bool check_overlap = overlap.Matrix().isApprox(overlap_ref, 0.0001);
  if (!check_overlap) {
    std::cout << std::endl;
    std::cout << "Ref" << std::endl;
    std::cout << overlap_ref << std::endl;
    std::cout << "Result" << std::endl;
    std::cout << overlap.Matrix() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_overlap, 1);
}

BOOST_AUTO_TEST_CASE(aocoulomb_inv_test) {
  QMMolecule mol = Methane();
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);

  AOCoulomb cou;
  cou.Fill(aobasis);

  Eigen::MatrixXd PseudoInvSqrt = cou.Pseudo_InvSqrt(1e-7);

  Eigen::MatrixXd Reformed = PseudoInvSqrt * PseudoInvSqrt * cou.Matrix();

  bool check_inv = Reformed.isApprox(Eigen::MatrixXd::Identity(17, 17), 0.0001);
  if (!check_inv) {
    std::cout << "reformed" << endl;
    std::cout << Reformed << endl;
  }
  BOOST_CHECK_EQUAL(check_inv, 1);
}

BOOST_AUTO_TEST_CASE(large_l_test) {

  QMMolecule mol("C", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/C2.xyz");

  BasisSet basisset;
  basisset.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/G.xml");

  BasisSet auxbasisset;
  auxbasisset.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/I.xml");
  AOBasis dftbasis;
  dftbasis.Fill(basisset, mol);

  AOBasis auxbasis;
  auxbasis.Fill(auxbasisset, mol);

  AOOverlap overlap;
  overlap.Fill(auxbasis);

  Eigen::MatrixXd overlap_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/overlap_ref_gi.mm");

  bool check_overlap = overlap.Matrix().isApprox(overlap_ref, 0.00001);

  BOOST_CHECK_EQUAL(check_overlap, 1);
  if (!check_overlap) {
    cout << "ref" << endl;
    cout << overlap_ref << endl;
    cout << "result" << endl;
    cout << overlap.Matrix() << endl;
  }

  AOCoulomb coulomb;
  coulomb.Fill(auxbasis);
  Eigen::MatrixXd coulomb_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/coulomb_ref_gi.mm");

  bool check_coulomb = coulomb.Matrix().isApprox(coulomb_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_coulomb, 1);
  if (!check_coulomb) {
    cout << "ref" << endl;
    cout << coulomb_ref << endl;
    cout << "result" << endl;
    cout << coulomb.Matrix() << endl;
  }

  AOKinetic kinetic;
  kinetic.Fill(dftbasis);

  Eigen::MatrixXd kinetic_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/kinetic_ref_gi.mm");

  bool check_kinetic = kinetic.Matrix().isApprox(kinetic_ref, 0.00001);

  BOOST_CHECK_EQUAL(check_kinetic, 1);
  if (!check_kinetic) {
    cout << "ref" << endl;
    cout << kinetic_ref << endl;
    cout << "result" << endl;
    cout << kinetic.Matrix() << endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
