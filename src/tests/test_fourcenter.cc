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

#define BOOST_TEST_MODULE fourcenter_test

// Third party includes
#include <boost/test/unit_test.hpp>

#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/fourcenter.h"
#include "votca/xtp/qmmolecule.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(fourcenter_test)

BOOST_AUTO_TEST_CASE(small_l_test) {

  QMMolecule mol(" ", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/fourcenter/molecule.xyz");

  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/fourcenter/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);

  const AOShell& shell0 = aobasis.getShell(2);
  const AOShell& shell1 = aobasis.getShell(3);
  const AOShell& shell2 = aobasis.getShell(4);
  const AOShell& shell4 = aobasis.getShell(5);
  FCMatrix fcenter;
  Eigen::Tensor<double, 4> block(shell0.getNumFunc(), shell4.getNumFunc(),
                                 shell2.getNumFunc(), shell1.getNumFunc());
  block.setZero();
  fcenter.FillFourCenterRepBlock(block, shell0, shell4, shell2, shell1);

  Eigen::Map<Eigen::VectorXd> mapped_result(block.data(), block.size());
  Eigen::VectorXd ref = Eigen::VectorXd::Zero(block.size());
  ref << 0.021578, 0.0112696, 0.0112696, 0.0112696, 0.021578, 0.0112696,
      0.0112696, 0.0112696, 0.021578;
  Eigen::TensorMap<Eigen::Tensor<double, 4> > ref_block(
      ref.data(), shell0.getNumFunc(), shell4.getNumFunc(), shell2.getNumFunc(),
      shell1.getNumFunc());

  bool check = mapped_result.isApprox(ref, 0.0001);
  BOOST_CHECK_EQUAL(check, 1);
  if (!check) {
    cout << "ref" << endl;
    cout << ref_block << endl;
    cout << "result" << endl;
    cout << block << endl;
  }
}

BOOST_AUTO_TEST_CASE(large_l_test) {

  QMMolecule mol("C", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/fourcenter/C2.xyz");

  BasisSet basisset;
  basisset.Load(std::string(XTP_TEST_DATA_FOLDER) + "/fourcenter/G.xml");

  AOBasis dftbasis;
  dftbasis.Fill(basisset, mol);

  FCMatrix fcenter;
  Eigen::Tensor<double, 4> block(
      dftbasis.getShell(0).getNumFunc(), dftbasis.getShell(1).getNumFunc(),
      dftbasis.getShell(0).getNumFunc(), dftbasis.getShell(1).getNumFunc());
  block.setZero();
  fcenter.FillFourCenterRepBlock(block, dftbasis.getShell(0),
                                 dftbasis.getShell(1), dftbasis.getShell(0),
                                 dftbasis.getShell(1));
  // we only check the first and last 600 values because this gets silly quite
  // quickly
  Eigen::Map<Eigen::VectorXd> mapped_result(block.data(), block.size());
  Eigen::VectorXd ref_head = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/fourcenter/largeLbegin.mm");

  Eigen::VectorXd ref_tail = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/fourcenter/largeLend.mm");

  bool check_head = mapped_result.head<600>().isApprox(ref_head, 0.0001);
  BOOST_CHECK_EQUAL(check_head, 1);
  if (!check_head) {
    cout << "ref" << endl;
    cout << ref_head.transpose() << endl;
    cout << "result" << endl;
    cout << mapped_result.head<600>().transpose() << endl;
  }

  bool check_tail = mapped_result.tail<600>().isApprox(ref_tail, 0.0001);
  BOOST_CHECK_EQUAL(check_tail, 1);
  if (!check_tail) {
    cout << "ref" << endl;
    cout << ref_tail.transpose() << endl;
    cout << "result" << endl;
    cout << mapped_result.tail<600>().transpose() << endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
