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

#define BOOST_TEST_MODULE aotransform_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <votca/xtp/aotransform.h>
#include <votca/xtp/orbitals.h>

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(aotransform_test)

BOOST_AUTO_TEST_CASE(transform) {
  QMAtom a(0, "C", Eigen::Vector3d::Zero());
  QMMolecule mol("zero", 0);
  mol.push_back(a);

  ofstream basisfile("all.xml");
  basisfile << "<basis name=\"test\">" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"P\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"D\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"D\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"F\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"F\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"G\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"G\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"H\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"H\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"I\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"I\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  BasisSet bs;
  bs.Load("all.xml");
  AOBasis basis;
  basis.Fill(bs, mol);

  std::array<Eigen::MatrixXd, 7> ref;
  ref[0] = Eigen::MatrixXd::Zero(1, 1);
  ref[1] = Eigen::MatrixXd::Zero(3, 3);
  ref[2] = Eigen::MatrixXd::Zero(6, 5);
  ref[3] = Eigen::MatrixXd::Zero(10, 7);
  ref[4] = Eigen::MatrixXd::Zero(15, 9);
  ref[5] = Eigen::MatrixXd::Zero(21, 11);
  ref[6] = Eigen::MatrixXd::Zero(28, 13);

  int ref_index = 0;
  for (const AOShell& shell : basis) {
    for (const AOGaussianPrimitive& gauss : shell) {
      Eigen::MatrixXd transform = AOTransform::getTrafo(gauss);
      bool check_transform = ref[ref_index].isApprox(transform, 1e-5);
      BOOST_CHECK_EQUAL(check_transform, 1);
      if (!check_transform) {
        std::cout << "ref " << shell.getType() << std::endl;
        std::cout << ref[ref_index] << std::endl;
        std::cout << "result" << std::endl;
        std::cout << transform << std::endl;
      }
      ref_index++;
    }
  }
}

BOOST_AUTO_TEST_CASE(composite_shell_transform) {
  QMAtom a(0, "C", Eigen::Vector3d::Zero());
  QMMolecule mol("zero", 0);
  mol.push_back(a);

  ofstream basisfile("composite_PD.xml");
  basisfile << "<basis name=\"test\">" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"PD\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"P\"/>"
            << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"D\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  BasisSet bs;
  bs.Load("composite_PD.xml");
  AOBasis basis;
  basis.Fill(bs, mol);

  Eigen::MatrixXd p = Eigen::MatrixXd::Zero(3, 3);

  Eigen::MatrixXd d = Eigen::MatrixXd::Zero(6, 5);

  Eigen::MatrixXd ref = Eigen::MatrixXd::Zero(9, 8);

  ref.topLeftCorner<3, 3>() = p;
  ref.bottomRightCorner<6, 5>() = d;

  for (const AOShell& shell : basis) {
    for (const AOGaussianPrimitive& gauss : shell) {
      Eigen::MatrixXd transform = AOTransform::getTrafo(gauss);
      bool check_transform = ref.isApprox(transform, 1e-5);
      BOOST_CHECK_EQUAL(check_transform, 1);
      if (!check_transform) {
        std::cout << "ref " << shell.getType() << std::endl;
        std::cout << ref << std::endl;
        std::cout << "result" << std::endl;
        std::cout << transform << std::endl;
      }
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
