/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE ewaldrealspaceinteractor_test

// Standard includes
#include <cmath>
#include <iostream>

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/ewaldrealspaceinteractor.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(ewaldrealspaceinteractor_test)

// A very small alpha makes erfc(alpha*r) ~= 1 for any r used in these
// tests, so it recovers the unscreened (bare Coulomb / bare multipole)
// expressions. This is the natural cross-check for the screened tensor:
// take away the screening almost entirely and it must fall back to
// ordinary electrostatics. 1e-13 (rather than a more modest value like
// 1e-8) is used because the leading residual is O(alpha*r) -- with the
// tight (1e-6%) tolerances used throughout this file, 1e-8 alone leaves
// a residual right at that tolerance's edge and was observed to fail.
constexpr double kTinyAlpha = 1e-13;

BOOST_AUTO_TEST_CASE(charge_charge_bare_coulomb_limit) {
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setCharge(1.0);
  StaticSite two(2, "H");
  two.setPos(Eigen::Vector3d::UnitX());
  two.setCharge(1.0);

  EwaldRealSpaceInteractor interactor(kTinyAlpha);
  double e = interactor.CalcStaticEnergy(one, two);

  BOOST_CHECK_CLOSE(e, 1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(charge_charge_explicit_erfc_formula) {
  const double alpha = 0.3;
  const double r = 2.5;

  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setCharge(1.5);
  StaticSite two(2, "H");
  two.setPos(r * Eigen::Vector3d::UnitZ());
  two.setCharge(-0.8);

  EwaldRealSpaceInteractor interactor(alpha);
  double e = interactor.CalcStaticEnergy(one, two);

  double e_ref = 1.5 * (-0.8) * std::erfc(alpha * r) / r;
  BOOST_CHECK_CLOSE(e, e_ref, 1e-9);
}

// Independent hand-derived reference for the screened charge-dipole /
// dipole-dipole energy, built directly from B0/B1/B2 rather than from
// EwaldRealSpaceInteractor itself, so this is a genuine cross-check and
// not a tautology.
double ScreenedMultipoleEnergyRef(double alpha, double q1,
                                  const Eigen::Vector3d& mu1, double q2,
                                  const Eigen::Vector3d& mu2,
                                  const Eigen::Vector3d& r_vec) {
  const double r = r_vec.norm();
  const double a = alpha;
  const double rSqrtPi = 0.5641895835477563;
  const double expTerm = rSqrtPi * std::exp(-a * a * r * r);
  const double B0 = std::erfc(a * r) / r;
  const double B1 = (B0 + 2.0 * a * expTerm) / (r * r);
  const double B2 = (3.0 * B1 + 4.0 * a * a * a * expTerm) / (r * r);

  const double mu1_r = mu1.dot(r_vec);
  const double mu2_r = mu2.dot(r_vec);
  const double phi1 = q1 * B0 + mu1_r * B1;
  // B2 already carries its (2l-1)=3 recursion factor, so this is
  // mu1_r*B2, not 3*mu1_r*B2 -- see EwaldRealSpaceInteractor::
  // EvaluateSource's own comment on the identical point.
  const Eigen::Vector3d E1 = q1 * B1 * r_vec + mu1_r * B2 * r_vec - B1 * mu1;
  (void)mu2_r;
  return q2 * phi1 - mu2.dot(E1);
}

BOOST_AUTO_TEST_CASE(charge_dipole_screened_reference) {
  const double alpha = 0.4;
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setCharge(0.6);
  StaticSite two(2, "H");
  two.setPos(Eigen::Vector3d(1.2, -0.3, 0.7));
  Vector9d mpoles2;
  mpoles2 << -0.6, 0.2, -0.5, 0.9, 0, 0, 0, 0, 0;
  two.setMultipole(mpoles2, 1);

  EwaldRealSpaceInteractor interactor(alpha);
  double e = interactor.CalcStaticEnergy(one, two);

  Eigen::Vector3d r_vec = two.getPos() - one.getPos();
  double e_ref = ScreenedMultipoleEnergyRef(
      alpha, one.getCharge(), Eigen::Vector3d::Zero(), two.getCharge(),
      two.getStaticDipole(), r_vec);

  BOOST_CHECK_CLOSE(e, e_ref, 1e-9);
}

BOOST_AUTO_TEST_CASE(dipole_dipole_unscreened_matches_textbook_formula) {
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  Vector9d mpoles1;
  mpoles1 << 0, 1, 1, 1, 0, 0, 0, 0, 0;
  one.setMultipole(mpoles1, 1);

  StaticSite two(2, "H");
  two.setPos(1.5 * Eigen::Vector3d::Ones());
  Vector9d mpoles2;
  mpoles2 << 0, -1, 1, 1, 0, 0, 0, 0, 0;
  two.setMultipole(mpoles2, 1);

  EwaldRealSpaceInteractor interactor(kTinyAlpha);
  double e = interactor.CalcStaticEnergy(one, two);

  // Standard unscreened dipole-dipole energy formula (same reference used
  // by test_eeinteractor.cc's own dipole test):
  // U = [ mu1.mu2 * r^2 - 3*(mu1.r)*(mu2.r) ] / r^5
  Eigen::Vector3d R = two.getPos() - one.getPos();
  Eigen::Vector3d mu1 = one.getStaticDipole();
  Eigen::Vector3d mu2 = two.getStaticDipole();
  double r = R.norm();
  double e_ref =
      (r * r * mu1.dot(mu2) - 3 * mu1.dot(R) * mu2.dot(R)) / std::pow(r, 5);

  BOOST_CHECK_CLOSE(e, e_ref, 1e-6);
}

BOOST_AUTO_TEST_CASE(apply_static_field_matches_energy_and_writes_v) {
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setCharge(2.0);

  PolarSite two(2, "H");
  two.setPos(2.0 * Eigen::Vector3d::UnitZ());

  EwaldRealSpaceInteractor interactor(0.35);
  double e = interactor.ApplyStaticField<StaticSite, Estatic::V>(one, two);

  // Field of a screened point charge along its own displacement axis:
  // E = q*B1*r_vec, purely along z here.
  const double r = 2.0;
  const double a = 0.35;
  const double rSqrtPi = 0.5641895835477563;
  const double expTerm = rSqrtPi * std::exp(-a * a * r * r);
  const double B0 = std::erfc(a * r) / r;
  const double B1 = (B0 + 2.0 * a * expTerm) / (r * r);
  Eigen::Vector3d field_ref = 2.0 * B1 * (r * Eigen::Vector3d::UnitZ());

  BOOST_CHECK(field_ref.isApprox(two.V(), 1e-9));
  // two has zero charge and zero static dipole, so the returned energy
  // (q2*phi1 - mu2.E1) must be exactly zero regardless of the field.
  BOOST_CHECK_CLOSE(e, 0.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(induced_thole_damping_short_vs_long_range) {
  // At tiny alpha, erfc screening is inactive, so any difference from the
  // plain (undamped) dipole-dipole energy is attributable entirely to
  // Thole damping.
  EwaldRealSpaceInteractor interactor(kTinyAlpha);

  PolarSite close1(1, "H");
  close1.setpolarization(Eigen::Matrix3d::Identity());
  close1.setPos(Eigen::Vector3d::Zero());
  close1.setInduced_Dipole(Eigen::Vector3d(0.3, 0, 0));

  PolarSite close2(2, "H");
  close2.setpolarization(Eigen::Matrix3d::Identity());
  close2.setPos(0.5 * Eigen::Vector3d::UnitX());
  close2.setInduced_Dipole(Eigen::Vector3d(0.2, 0, 0));

  double e_close = interactor.CalcInducedEnergy(close1, close2);

  Eigen::Vector3d R = close2.getPos() - close1.getPos();
  double r = R.norm();
  Eigen::Vector3d mu1 = close1.getInducedDipole();
  Eigen::Vector3d mu2 = close2.getInducedDipole();
  double e_undamped_close =
      (r * r * mu1.dot(mu2) - 3 * mu1.dot(R) * mu2.dot(R)) / std::pow(r, 5);

  // Short range (au3 well under 40 for unit polarizability and this
  // separation): Thole damping must be active, i.e. the screened energy
  // must differ measurably from the undamped point-dipole formula.
  BOOST_CHECK(std::abs(e_close - e_undamped_close) > 1e-3);

  PolarSite far1(3, "H");
  far1.setpolarization(Eigen::Matrix3d::Identity());
  far1.setPos(Eigen::Vector3d::Zero());
  far1.setInduced_Dipole(Eigen::Vector3d(0.3, 0, 0));

  PolarSite far2(4, "H");
  far2.setpolarization(Eigen::Matrix3d::Identity());
  far2.setPos(20.0 * Eigen::Vector3d::UnitX());
  far2.setInduced_Dipole(Eigen::Vector3d(0.2, 0, 0));

  double e_far = interactor.CalcInducedEnergy(far1, far2);

  Eigen::Vector3d Rf = far2.getPos() - far1.getPos();
  double rf = Rf.norm();
  Eigen::Vector3d mu1f = far1.getInducedDipole();
  Eigen::Vector3d mu2f = far2.getInducedDipole();
  double e_undamped_far = (rf * rf * mu1f.dot(mu2f) -
                           3 * mu1f.dot(Rf) * mu2f.dot(Rf)) /
                          std::pow(rf, 5);

  // Long range: Thole damping factors must have saturated to 1, so the
  // screened (tiny alpha) energy must match the undamped point-dipole
  // formula closely.
  BOOST_CHECK_CLOSE(e_far, e_undamped_far, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
