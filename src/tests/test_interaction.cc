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
  Eigen::VectorXd multipole = Eigen::VectorXd::Zero(9);
  multipole << 1, 2, 3, 4, 8, 7, 2, 3.3, -0.5;
  ps.setMultipole(multipole);
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

  Eigen::VectorXd multipoles = Eigen::VectorXd::Zero(9);
  multipoles << 1, 1, 0, 0, 0, 1, 0, 0,
      0;  // q=1, mu_x=1 and Q_21c=1 the rest is 0
  ps.setMultipole(multipoles);
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

BOOST_AUTO_TEST_CASE(interaction_test) {
  PolarSite ps1(1, "ps1");
  PolarSite ps2(2, "ps2", Eigen::Vector3d::UnitX());

  Eigen::VectorXd mp1 = Eigen::VectorXd::Zero(1);
  Eigen::VectorXd mp2 = Eigen::VectorXd::Zero(1);
  mp1 << 1;
  mp2 << -1;
  ps1.setPolarisable(false);
  ps2.setPolarisable(false);
  ps1.setMultipole(mp1);
  ps2.setMultipole(mp2);

  double Energyref = -1;
  double Energy = ps1.InteractStatic(ps2);
  BOOST_CHECK_EQUAL(std::abs(Energy - Energyref) < 1e-9, true);

  bool check_field = ps1.getField().isApprox(ps2.getField(), 1e-5);
  if (!check_field) {
    std::cout << "Field at ps1" << std::endl;
    std::cout << ps1.getField() << std::endl;
    std::cout << "Field at ps2" << std::endl;
    std::cout << ps2.getField() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_field, true);

  bool check_potential =
      std::abs(ps1.getPotential() + ps2.getPotential()) < 1e-5;
  if (!check_potential) {
    std::cout << "Potential at ps1" << std::endl;
    std::cout << ps1.getPotential() << std::endl;
    std::cout << "Potential at ps2" << std::endl;
    std::cout << ps2.getPotential() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_potential, true);
  PolarSite ps3(3, "ps3");
  PolarSite ps4(4, "ps4", Eigen::Vector3d::UnitZ());
  Eigen::VectorXd multipole = Eigen::VectorXd::Zero(9);
  multipole << 1, 2, 3, 4, 8, 7, 2, 3.3, -0.5;

  ps3.setPolarisable(true);
  ps4.setPolarisable(true);
  ps3.setMultipole(multipole);
  ps4.setMultipole(multipole);
  ps3.InteractStatic(ps4);
}

BOOST_AUTO_TEST_CASE(induction_test) {
  PolarSite ps1(1, "ps1");
  PolarSite ps2(2, "ps2", Eigen::Vector3d::UnitX());

  Eigen::VectorXd mp1 = Eigen::VectorXd::Zero(1);
  Eigen::VectorXd mp2 = Eigen::VectorXd::Zero(1);
  mp1 << 1;
  mp2 << -1;
  ps1.setPolarisable(true);
  ps2.setPolarisable(true);
  ps1.setMultipole(mp1);
  ps2.setMultipole(mp2);
  Eigen::Matrix3d poltensor = Eigen::Matrix3d::Zero();
  poltensor << 2, 1, 0, 1, 3, 1, 0, 1, 2.5;

  ps1.setPolarisation(poltensor);
  ps2.setPolarisation(poltensor);

  // double Energy= ps1.InteractStatic(ps2);
  ps1.Induce(1);
  ps2.Induce(1);
  // double alpha=0.39;
  // double InductionEnergy=ps1.InteractInduction(ps2,alpha);
}

BOOST_AUTO_TEST_SUITE_END()
