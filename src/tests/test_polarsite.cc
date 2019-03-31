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

#define BOOST_TEST_MODULE polararsite_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/polarsite.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(polararsite_test)

BOOST_AUTO_TEST_CASE(constructors_test) { PolarSite ps(1, "ps1"); }

BOOST_AUTO_TEST_CASE(getters_test) {
  PolarSite ps(1, "ps2");
  BOOST_CHECK_EQUAL(ps.getId(), 1);
  BOOST_CHECK_EQUAL(ps.getElement(), "ps2");
}

BOOST_AUTO_TEST_CASE(multipole_test) {
  PolarSite ps(1, "ps2");
  Vector9d multipole = Vector9d::Zero(9);
  multipole << 1, 2, 3, 4, 8, 7, 2, 3.3, -0.5;
  ps.setMultipole(multipole, 2);
  bool check_mpoles = multipole.isApprox(ps.getPermMultipole(), 0.0001);
  BOOST_CHECK_EQUAL(check_mpoles, true);

  bool check_rank = (ps.getRank() == 2);
  BOOST_CHECK_EQUAL(check_rank, true);
}

BOOST_AUTO_TEST_CASE(translate_test) {
  PolarSite ps(1, "ps2");
  Eigen::Vector3d shift;
  shift << 0, 0, 5;
  ps.Translate(shift);
  BOOST_CHECK_EQUAL(shift.isApprox(ps.getPos(), 1e-5), true);
}

BOOST_AUTO_TEST_CASE(rotation_test) {
  PolarSite ps(1, "ps2", Eigen::Vector3d::UnitY());

  Eigen::Matrix3d R = Eigen::Matrix3d::Zero();  // Rotation around z axes
  R << 0, -1, 0, 1, 0, 0, 0, 0, 1;

  Vector9d multipoles = Vector9d::Zero(9);
  multipoles << 1, 1, 0, 0, 0, 1, 0, 0,
      0;  // q=1, mu_x=1 and Q_21c=1 the rest is 0
  ps.setMultipole(multipoles, 2);
  ps.Rotate(R, Eigen::Vector3d::Zero());
  bool equalpos = ps.getPos().isApprox(Eigen::Vector3d(-1, 0, 0), 1e-5);
  if (!equalpos) {
    std::cout << "Result " << std::endl;
    std::cout << ps.getPos() << std::endl;
    std::cout << "Reference" << std::endl;
    std::cout << Eigen::Vector3d(-1, 0, 0) << std::endl;
  }
  BOOST_CHECK_EQUAL(equalpos, true);

  Eigen::VectorXd rotmultipoles = Eigen::VectorXd::Zero(9);
  rotmultipoles << 1, 0, 1, 0, 0, 0, 1, 0, 0;  // q=1, mu_y=1 and Q_21s=1 is 0
  bool equalmultipoles = rotmultipoles.isApprox(ps.getPermMultipole(), 1e-5);
  if (!equalmultipoles) {
    std::cout << "Result " << std::endl;
    std::cout << ps.getPermMultipole() << std::endl;
    std::cout << "Reference" << std::endl;
    std::cout << rotmultipoles << std::endl;
  }
  BOOST_CHECK_EQUAL(equalmultipoles, true);
}

BOOST_AUTO_TEST_SUITE_END()
