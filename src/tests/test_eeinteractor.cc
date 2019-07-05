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
#include <votca/xtp/dipoledipoleinteraction.h>
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

BOOST_AUTO_TEST_CASE(static_case_monopole_field) {

  Vector9d mpoles1;
  mpoles1 << -3, 0, 0, 0, 0, 0, 0, 0, 0;

  StaticSegment seg1("one", 1);
  PolarSegment seg2("two", 2);
  StaticSite one(1, "H");
  one.setMultipole(mpoles1, 0);
  one.setPos(Eigen::Vector3d::Zero());

  PolarSite two(2, "H");
  two.setPos(3 * Eigen::Vector3d::UnitZ());

  seg1.push_back(one);
  seg2.push_back(two);

  eeInteractor interactor;
  interactor.ApplyStaticField<StaticSegment, Estatic::V>(seg1, seg2);

  Eigen::Vector3d field_ref;
  field_ref << 0, 0, 0.33333333;
  bool field_check = field_ref.isApprox(seg2[0].V(), 1e-6);
  BOOST_CHECK_EQUAL(field_check, true);
  if (!field_check) {
    std::cout << "ref" << std::endl;
    std::cout << field_ref << std::endl;
    std::cout << "field" << std::endl;
    std::cout << seg2[0].V() << std::endl;
  }

  interactor.ApplyStaticField<StaticSegment, Estatic::noE_V>(seg1, seg2);

  bool field_check_2 = field_ref.isApprox(seg2[0].V_noE(), 1e-6);
  BOOST_CHECK_EQUAL(field_check, true);
  if (!field_check_2) {
    std::cout << "ref" << std::endl;
    std::cout << field_ref << std::endl;
    std::cout << "field" << std::endl;
    std::cout << seg2[0].V_noE() << std::endl;
  }
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

BOOST_AUTO_TEST_CASE(static_case_full_tensor_1) {
  eeInteractor interactor;
  Eigen::Matrix<double, 9, 9> result;
  for (int i = 0; i < 9; i++) {
    Vector9d mpoles1 = Vector9d::Zero();
    mpoles1(i) = 1.0;
    StaticSegment seg1("one", 1);
    StaticSite one(1, "H");
    one.setPos(Eigen::Vector3d::Zero());
    one.setMultipole(mpoles1, 2);
    seg1.push_back(one);
    for (int j = 0; j < 9; j++) {
      Vector9d mpoles2 = Vector9d::Zero();
      mpoles2(j) = 1.0;
      StaticSegment seg2("two", 2);
      StaticSite two(2, "H");
      two.setPos(3 * Eigen::Vector3d::Ones());
      two.setMultipole(mpoles2, 2);
      seg2.push_back(two);
      result(i, j) = interactor.CalcStaticEnergy(seg1, seg2);
    }
  }

  Eigen::Matrix<double, 9, 9> ref = Eigen::Matrix<double, 9, 9>::Zero();
  ref << 0.19245, -0.0213833, -0.0213833, -0.0213833, 0, 0.00411523, 0.00411523,
      0, 0.00411523, 0.0213833, 8.67362e-19, -0.00712778, -0.00712778,
      0.000791976, 0.000914495, 0.00228624, -0.00137174, 0.000914495, 0.0213833,
      -0.00712778, 8.67362e-19, -0.00712778, 0.000791976, 0.00228624,
      0.000914495, 0.00137174, 0.000914495, 0.0213833, -0.00712778, -0.00712778,
      8.67362e-19, -0.00158395, 0.000914495, 0.000914495, 0, 0.00228624, 0,
      -0.000791976, -0.000791976, 0.00158395, -0.000615981, -0.000254026,
      -0.000254026, 0, 0.000508053, 0.00411523, -0.000914495, -0.00228624,
      -0.000914495, -0.000254026, 0.000410654, 0.000586649, -0.000439986,
      0.000586649, 0.00411523, -0.00228624, -0.000914495, -0.000914495,
      -0.000254026, 0.000586649, 0.000410654, 0.000439986, 0.000586649, 0,
      0.00137174, -0.00137174, 0, 0, -0.000439986, 0.000439986, -0.000615981, 0,
      0.00411523, -0.000914495, -0.000914495, -0.00228624, 0.000508053,
      0.000586649, 0.000586649, 0, 0.000410654;

  bool tensor_check = result.isApprox(ref, 1e-6);
  BOOST_CHECK_EQUAL(tensor_check, true);
  if (!tensor_check) {
    std::cout << "result" << std::endl;
    std::cout << result << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(static_case_full_tensor_2) {
  eeInteractor interactor;
  Eigen::Matrix<double, 9, 9> result;
  for (int i = 0; i < 9; i++) {
    Vector9d mpoles1 = Vector9d::Zero();
    mpoles1(i) = 1.0;
    StaticSegment seg1("one", 1);
    StaticSite one(1, "H");
    one.setPos(Eigen::Vector3d::Zero());
    one.setMultipole(mpoles1, 2);
    seg1.push_back(one);
    for (int j = 0; j < 9; j++) {
      Vector9d mpoles2 = Vector9d::Zero();
      mpoles2(j) = 1.0;
      StaticSegment seg2("two", 2);
      StaticSite two(2, "H");
      two.setPos(-2.1 * Eigen::Vector3d::Ones() + Eigen::Vector3d::UnitX());
      two.setMultipole(mpoles2, 2);
      seg2.push_back(two);
      result(i, j) = interactor.CalcStaticEnergy(seg1, seg2);
    }
  }

  Eigen::Matrix<double, 9, 9> ref = Eigen::Matrix<double, 9, 9>::Zero();
  ref << 0.315754, 0.0346291, 0.0661101, 0.0661101, 0.0050219, 0.012558,
      0.0239744, -0.00869818, 0.012558, -0.0346291, 0.0200876, -0.0217511,
      -0.0217511, -0.00620633, 0.00453012, -0.0131465, 0.0107497, 0.00453012,
      -0.0661101, -0.0217511, -0.0100438, -0.0415248, -0.0118485, -0.0131465,
      -0.0136814, -0.00231059, -0.00716646, -0.0661101, -0.0217511, -0.0415248,
      -0.0100438, 0.00792526, -0.00716646, -0.0136814, 0.00910577, -0.0131465,
      0.0050219, 0.00620633, 0.0118485, -0.00792526, -0.00806073, 0.000243418,
      0.000464707, -0.00450468, 0.00650363, 0.012558, -0.00453012, 0.0131465,
      0.00716646, 0.000243418, 0.000172264, 0.00750975, -0.00765029, 0.0143368,
      0.0239744, 0.0131465, 0.0136814, 0.0136814, 0.000464707, 0.00750975,
      0.0105754, -0.000804897, 0.00750975, -0.00869818, -0.0107497, 0.00231059,
      -0.00910577, -0.00450468, -0.00765029, -0.000804897, -0.00285918,
      -0.00403595, 0.012558, -0.00453012, 0.00716646, 0.0131465, 0.00650363,
      0.0143368, 0.00750975, -0.00403595, 0.000172264;

  bool tensor_check = result.isApprox(ref, 1e-5);
  BOOST_CHECK_EQUAL(tensor_check, true);
  if (!tensor_check) {
    std::cout << "result" << std::endl;
    std::cout << result << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(polar_case_monopole) {

  Vector9d mpoles1;
  mpoles1 << 1, 0, 0, 0, 0, 0, 0, 0, 0;

  PolarSegment seg1("one", 1);
  PolarSegment seg2("two", 2);
  PolarSite one(1, "H");
  one.setPolarisation(Eigen::Matrix3d::Identity());
  one.setMultipole(mpoles1, 0);
  one.setPos(Eigen::Vector3d::Zero());

  PolarSite two(2, "H");
  two.setPolarisation(Eigen::Matrix3d::Identity());
  two.setMultipole(mpoles1, 0);
  two.setPos(Eigen::Vector3d::UnitZ());

  seg1.push_back(one);
  seg2.push_back(two);
  double exp_damp =
      std::numeric_limits<double>::max();  // should disable thole damping
  eeInteractor interactor(exp_damp);
  double estat =
      interactor.ApplyStaticField<PolarSegment, Estatic::noE_V>(seg1, seg2);
  interactor.ApplyStaticField<PolarSegment, Estatic::noE_V>(seg2, seg1);
  BOOST_CHECK_CLOSE(estat, 1, 1e-12);
  Eigen::Vector3d field_ref;
  field_ref << 0, 0, -1;
  bool field_check = field_ref.isApprox(seg2[0].V_noE(), 1e-6);
  BOOST_CHECK_EQUAL(field_check, true);
  if (!field_check) {
    std::cout << "ref" << std::endl;
    std::cout << field_ref << std::endl;
    std::cout << "field" << std::endl;
    std::cout << seg2[0].V_noE() << std::endl;
  }

  bool field_check_2 = field_ref.isApprox(-seg1[0].V_noE(), 1e-6);
  BOOST_CHECK_EQUAL(field_check, true);
  if (!field_check_2) {
    std::cout << "ref" << std::endl;
    std::cout << -field_ref << std::endl;
    std::cout << "field" << std::endl;
    std::cout << seg1[0].V_noE() << std::endl;
  }

  std::vector<PolarSegment> segments = {seg1, seg2};
  Eigen::VectorXd b = Eigen::VectorXd::Zero(6);
  int index = 0;
  for (PolarSegment& seg : segments) {
    for (const PolarSite& site : seg) {
      const Eigen::Vector3d V = site.V() + site.V_noE();
      b.segment<3>(index) = -V;
      index += 3;
    }
  }

  DipoleDipoleInteraction A(interactor, segments);
  Eigen::ConjugateGradient<DipoleDipoleInteraction, Eigen::Lower | Eigen::Upper>
      cg;
  cg.setMaxIterations(100);
  cg.setTolerance(1e-9);
  cg.compute(A);
  Eigen::VectorXd x = cg.solve(b);
  index = 0;
  for (PolarSegment& seg : segments) {
    for (PolarSite& site : seg) {
      site.setInduced_Dipole(x.segment<3>(index));
      index += 3;
    }
  }

  Eigen::Vector3d dipole_ref;
  dipole_ref << 0, 0, 1.0 / 3.0;

  bool dipole_check =
      dipole_ref.isApprox(-segments[0][0].Induced_Dipole(), 1e-6);
  BOOST_CHECK_EQUAL(field_check, true);
  if (!dipole_check) {
    std::cout << "ref" << std::endl;
    std::cout << -dipole_ref << std::endl;
    std::cout << "dipole" << std::endl;
    std::cout << segments[0][0].Induced_Dipole() << std::endl;
  }

  bool dipole_check2 =
      dipole_ref.isApprox(segments[1][0].Induced_Dipole(), 1e-6);
  if (!dipole_check2) {
    std::cout << "ref" << std::endl;
    std::cout << dipole_ref << std::endl;
    std::cout << "dipole" << std::endl;
    std::cout << segments[1][0].Induced_Dipole() << std::endl;
  }

  eeInteractor::E_terms epolar =
      interactor.CalcPolarEnergy(segments[0], segments[1]);
  double einternal = segments[0][0].InternalEnergy();
  einternal += segments[1][0].InternalEnergy();

  BOOST_CHECK_CLOSE(epolar.E_indu_indu(), 2.0 / 9.0, 1e-12);
  BOOST_CHECK_CLOSE(epolar.E_indu_stat(), -2.0 / 3.0, 1e-12);
  BOOST_CHECK_CLOSE(einternal, 1.0 / 9.0, 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
