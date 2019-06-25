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

#define BOOST_TEST_MODULE eeinteractor_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/xtp/eeinteractor.h>
using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(eeinteractor_test)

BOOST_AUTO_TEST_CASE(static_case_charges) {

  StaticSegment seg1("one", 1);
  StaticSegment seg2("two", 2);
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setCharge(1);
  StaticSite two(2, "H");
  two.setPos(Eigen::Vector3d::UnitX());
  two.setCharge(1);

  seg1.push_back(one);
  seg2.push_back(two);

  eeInteractor interactor;
  double energy_field = interactor.CalcStaticEnergy(seg1, seg2);

  BOOST_CHECK_CLOSE(energy_field, 1, 1e-12);

  seg2[0].setPos(Eigen::Vector3d::UnitZ());
  double energy_field2 = interactor.CalcStaticEnergy(seg1, seg2);

  BOOST_CHECK_CLOSE(energy_field, energy_field2, 1e-12);
}

double CalcDipoleEnergy(StaticSegment& seg1, StaticSegment& seg2) {
  Eigen::Vector3d R = (seg2[0].getPos() - seg1[0].getPos());
  Eigen::Vector3d mu1 = seg1[0].getDipole();
  Eigen::Vector3d mu2 = seg2[0].getDipole();
  double r = R.norm();
  double e_ref =
      (r * r * mu1.dot(mu2) - 3 * mu1.dot(R) * mu2.dot(R)) / std::pow(r, 5);
  return e_ref;
}

BOOST_AUTO_TEST_CASE(static_case_dipoles) {

  Vector9d mpoles1;
  mpoles1 << 0, 1, 1, 1, 0, 0, 0, 0, 0;
  Vector9d mpoles2;
  mpoles2 << 0, -1, 1, 1, 0, 0, 0, 0, 0;

  StaticSegment seg1("one", 1);
  StaticSegment seg2("two", 2);
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setMultipole(mpoles1, 1);
  StaticSite two(2, "H");
  two.setPos(1.5 * Eigen::Vector3d::Ones());
  two.setMultipole(mpoles2, 1);

  seg1.push_back(one);
  seg2.push_back(two);

  eeInteractor interactor;
  double energy_field = interactor.CalcStaticEnergy(seg1, seg2);

  // explicit formula
  double e_ref = CalcDipoleEnergy(seg1, seg2);
  BOOST_CHECK_CLOSE(energy_field, e_ref, 1e-12);

  seg2[0].setPos(Eigen::Vector3d::UnitZ());
  double energy_field2 = interactor.CalcStaticEnergy(seg1, seg2);
  double e_ref2 = CalcDipoleEnergy(seg1, seg2);
  BOOST_CHECK_CLOSE(energy_field2, e_ref2, 1e-12);

  double energy_field3 = interactor.CalcStaticEnergy(seg2, seg1);
  double e_ref3 = CalcDipoleEnergy(seg2, seg1);
  BOOST_CHECK_CLOSE(e_ref2, e_ref3, 1e-12);
  BOOST_CHECK_CLOSE(energy_field3, e_ref3, 1e-12);
  BOOST_CHECK_CLOSE(energy_field3, energy_field2, 1e-12);
}

double CalcQuadrupoleEnergy(StaticSegment& seg1, StaticSegment& seg2) {
  Eigen::Vector3d R = (seg2[0].getPos() - seg1[0].getPos());
  double q20_1 = seg1[0].Q()[4];
  double q20_2 = seg2[0].Q()[4];
  double r = R.norm();
  Eigen::Vector3d Rhat = R / r;
  Eigen::Vector3d na = Eigen::Vector3d::UnitZ();
  Eigen::Vector3d nb = Eigen::Vector3d::UnitZ();
  double naR = Rhat.dot(na);
  double nbR = Rhat.dot(nb);
  // stone book p51 eqn 3.3.17
  double e_ref = q20_1 * q20_2 / std::pow(r, 5) * 0.75 *
                 (35 * naR * naR * nbR * nbR - 5 * naR - 5 * nbR -
                  20 * naR * nbR * na.dot(nb) + 2 * na.dot(nb) + 1);
  return e_ref;
}

BOOST_AUTO_TEST_CASE(static_case_quadrupoles) {

  Vector9d mpoles1;
  mpoles1 << 0, 0, 0, 0, 1, 0, 0, 0, 0;
  Vector9d mpoles2;
  mpoles2 << 0, 0, 0, 0, 1, 0, 0, 0, 0;

  StaticSegment seg1("one", 1);
  StaticSegment seg2("two", 2);
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setMultipole(mpoles1, 2);
  StaticSite two(2, "H");
  // for the quadrupole formulathe second molecule must be along the z axis
  two.setPos(3 * Eigen::Vector3d::UnitZ());
  two.setMultipole(mpoles2, 2);

  seg1.push_back(one);
  seg2.push_back(two);

  eeInteractor interactor;
  double energy_field = interactor.CalcStaticEnergy(seg1, seg2);

  // explicit formula
  double e_ref = CalcQuadrupoleEnergy(seg1, seg2);
  BOOST_CHECK_CLOSE(energy_field, e_ref, 1e-12);
}

BOOST_AUTO_TEST_CASE(polar_case_quadrupoles_field) {

  Vector9d mpoles1;
  mpoles1 << -1, 1, 0, 0, 1, 1, 1, 1, 1;

  StaticSegment seg1("one", 1);
  PolarSegment seg2("two", 2);
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());

  PolarSite two(2, "H");
  // for the quadrupole formulathe second molecule must be along the z axis
  two.setPos(3 * Eigen::Vector3d::UnitZ());

  seg1.push_back(one);
  seg2.push_back(two);

  eeInteractor interactor;
  interactor.ApplyStaticField(seg1, seg2);
}

BOOST_AUTO_TEST_CASE(static_case_quadrupoles_orientation) {
  // using Stone page 53 fig 3.3
  Vector9d mpoles1;
  mpoles1 << 0, 0, 0, 0, 1, 0, 0, 0, 0;
  Vector9d mpoles2;
  mpoles2 << 0, 0, 0, 0, 1, 0, 0, 0, 0;

  StaticSegment seg1("one", 1);
  StaticSegment seg2("two", 2);
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setMultipole(mpoles1, 2);
  StaticSite two(2, "H");
  // for the quadrupole formulathe second molecule must be along the z axis
  two.setPos(Eigen::Vector3d::UnitZ());
  two.setMultipole(mpoles2, 2);

  seg1.push_back(one);
  seg2.push_back(two);

  eeInteractor interactor;
  double energy_field = interactor.CalcStaticEnergy(seg1, seg2);
  // config a
  BOOST_CHECK_CLOSE(energy_field, 6, 1e-12);
  StaticSite& site1 = seg1[0];
  StaticSite& site2 = seg2[0];
  // site1 is at origin
  Eigen::Matrix3d rot;
  rot << 1, 0, 0, 0, 0, 1, 0, 1, 0;
  site1.Rotate(rot, Eigen::Vector3d::Zero());
  energy_field = interactor.CalcStaticEnergy(seg1, seg2);
  BOOST_CHECK_CLOSE(energy_field, -3, 1e-12);
  // config c
  site2.Rotate(rot, Eigen::Vector3d::UnitZ());
  energy_field = interactor.CalcStaticEnergy(seg1, seg2);
  BOOST_CHECK_CLOSE(energy_field, 2.25, 1e-12);
  // config d
  Eigen::Matrix3d rot2;
  rot2 << 0, 1, 0, 1, 0, 0, 0, 0, 1;
  site1.Rotate(rot2, Eigen::Vector3d::Zero());
  energy_field = interactor.CalcStaticEnergy(seg1, seg2);
  BOOST_CHECK_CLOSE(energy_field, 0.75, 1e-12);

  // reset to config a
  site1.setMultipole(mpoles1, 2);
  site2.setMultipole(mpoles2, 2);
  Eigen::Matrix3d rot3;
  double sqrt2 = 1 / std::sqrt(2);
  rot3 << 1, 0, 0, 0, sqrt2, sqrt2, 0, -sqrt2, sqrt2;
  site2.Rotate(rot3, Eigen::Vector3d::UnitZ());
  site1.Rotate(rot3, Eigen::Vector3d::Zero());
  energy_field = interactor.CalcStaticEnergy(seg1, seg2);
  double ref = -(2.0 + 7.0 / 16.0);
  BOOST_CHECK_CLOSE(energy_field, ref, 1e-12);
}

BOOST_AUTO_TEST_CASE(static_case_quadrupoles_dipoles_orientation) {
  // using Stone page 55 fig 3.5
  Vector9d mpoles1;
  mpoles1 << 0, 0, 0, 0, 1, 0, 0, 0, 0;
  Vector9d mpoles2;
  mpoles2 << 0, 0, 0, 1, 0, 0, 0, 0, 0;

  StaticSegment seg1("one", 1);
  StaticSegment seg2("two", 2);
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setMultipole(mpoles1, 2);
  StaticSite two(2, "H");
  // site 2 is above site 1
  two.setPos(Eigen::Vector3d::UnitZ());
  two.setMultipole(mpoles2, 2);

  seg1.push_back(one);
  seg2.push_back(two);

  eeInteractor interactor;
  double energy_field = interactor.CalcStaticEnergy(seg1, seg2);
  double energy_field_rev = interactor.CalcStaticEnergy(seg2, seg1);
  // config a
  BOOST_CHECK_CLOSE(energy_field, -3, 1e-12);
  BOOST_CHECK_CLOSE(energy_field, energy_field_rev, 1e-12);
  // config b
  StaticSite& site1 = seg1[0];
  StaticSite& site2 = seg2[0];
  // site1 is at origin
  Eigen::Matrix3d rot;
  rot << 1, 0, 0, 0, 1, 0, 0, 0, -1;
  site2.Rotate(rot, Eigen::Vector3d::UnitZ());
  energy_field = interactor.CalcStaticEnergy(seg1, seg2);
  energy_field_rev = interactor.CalcStaticEnergy(seg2, seg1);
  BOOST_CHECK_CLOSE(energy_field, 3, 1e-12);
  BOOST_CHECK_CLOSE(energy_field, energy_field_rev, 1e-12);
  // config c
  Eigen::Matrix3d rot2;
  rot2 << 1, 0, 0, 0, 0, 1, 0, 1, 0;
  site1.Rotate(rot2, Eigen::Vector3d::Zero());
  energy_field = interactor.CalcStaticEnergy(seg2, seg1);
  energy_field_rev = interactor.CalcStaticEnergy(seg2, seg1);
  BOOST_CHECK_CLOSE(energy_field, -1.5, 1e-12);
  BOOST_CHECK_CLOSE(energy_field, energy_field_rev, 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
