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
#include "votca/xtp/basisset.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE aotransform_test

// Third party includes
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aotransform.h"
#include "votca/xtp/orbitals.h"
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(aotransform_test)

BOOST_AUTO_TEST_CASE(transform_test) {
  QMAtom a(0, "C", Eigen::Vector3d::Zero());
  QMMolecule mol("zero", 0);
  mol.push_back(a);

  BasisSet bs;
  bs.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aotransform/all.xml");
  AOBasis basis;
  basis.Fill(bs, mol);

  std::array<Eigen::MatrixXd, 5> ref;

  for (unsigned i = 0; i < ref.size(); i++) {
    ref[i] = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/aotransform/ref_" +
        std::to_string(i) + ".mm");
  }

  Index ref_index = 0;
  for (const AOShell& shell : basis) {
    for (const AOGaussianPrimitive& gauss : shell) {
      Eigen::MatrixXd transform = AOTransform::getTrafo(gauss);
      bool check_transform = ref[ref_index].isApprox(transform, 1e-5);
      BOOST_CHECK_EQUAL(check_transform, 1);
      if (!check_transform) {
        std::cout << "ref " << xtp::EnumToString(shell.getL()) << std::endl;
        std::cout << ref[ref_index] << std::endl;
        std::cout << "result" << std::endl;
        std::cout << transform << std::endl;
      }
      ref_index++;
    }
  }
}

BOOST_AUTO_TEST_CASE(xintegrate) {

  BOOST_REQUIRE_THROW(AOTransform::XIntegrate(0, 0.1), std::runtime_error);

  BOOST_REQUIRE_THROW(AOTransform::XIntegrate(1, -0.1);, std::runtime_error);

  Eigen::VectorXd res1 = AOTransform::XIntegrate(1, 0.1);
  Eigen::VectorXd res1_ref = Eigen::VectorXd::Zero(1);
  res1_ref << 0.967643;
  bool check_res1 = res1.isApprox(res1_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_res1, 1);
  if (!check_res1) {
    std::cout << "ref" << std::endl;
    std::cout << res1_ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << res1 << std::endl;
  }
  Eigen::VectorXd res2 = AOTransform::XIntegrate(5, 0.1);
  Eigen::VectorXd res2_ref = Eigen::VectorXd::Zero(5);
  res2_ref << 0.967643, 0.314029, 0.186255, 0.132188, 0.102394;

  bool check_res2 = res2.isApprox(res2_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_res2, 1);
  if (!check_res1) {
    std::cout << "ref" << std::endl;
    std::cout << res2_ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << res2 << std::endl;
  }

  Eigen::VectorXd res3 = AOTransform::XIntegrate(1, 1e-12);
  Eigen::VectorXd res3_ref = Eigen::VectorXd::Zero(1);
  res3_ref[0] = 1.0;
  bool check_res3 = res3.isApprox(res3_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_res3, 1);
  if (!check_res3) {
    std::cout << "ref" << std::endl;
    std::cout << res3_ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << res3 << std::endl;
  }

  Eigen::VectorXd res4 = AOTransform::XIntegrate(5, 1e-12);
  Eigen::VectorXd res4_ref = Eigen::VectorXd::Zero(5);
  res4_ref << 1, 0.333333, 0.2, 0.142857, 0.111111;

  bool check_res4 = res4.isApprox(res4_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_res4, 1);
  if (!check_res4) {
    std::cout << "ref" << std::endl;
    std::cout << res4_ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << res4 << std::endl;
  }

  Eigen::VectorXd res5 = AOTransform::XIntegrate(1, 15);
  Eigen::VectorXd res5_ref = Eigen::VectorXd::Zero(1);
  res5_ref << 0.228823;
  bool check_res5 = res5.isApprox(res5_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_res5, 1);
  if (!check_res5) {
    std::cout << "ref" << std::endl;
    std::cout << res5_ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << res5 << std::endl;
  }

  Eigen::VectorXd res6 = AOTransform::XIntegrate(5, 15);
  Eigen::VectorXd res6_ref = Eigen::VectorXd::Zero(5);
  res6_ref << 0.228823, 0.00762742, 0.000762731, 0.000127112, 2.96492e-05;
  bool check_res6 = res6.isApprox(res6_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_res6, 1);
  if (!check_res6) {
    std::cout << "ref" << std::endl;
    std::cout << res6_ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << res6 << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(blocksize) {
  BOOST_CHECK_EQUAL(AOTransform::getBlockSize(0), 1);

  BOOST_CHECK_EQUAL(AOTransform::getBlockSize(0), 1);
  BOOST_CHECK_EQUAL(AOTransform::getBlockSize(1), 4);
  BOOST_CHECK_EQUAL(AOTransform::getBlockSize(2), 10);
  BOOST_CHECK_EQUAL(AOTransform::getBlockSize(3), 20);
}

BOOST_AUTO_TEST_SUITE_END()
