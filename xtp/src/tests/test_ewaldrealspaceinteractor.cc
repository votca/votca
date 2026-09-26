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
#include <vector>

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
  double e_undamped_far =
      (rf * rf * mu1f.dot(mu2f) - 3 * mu1f.dot(Rf) * mu2f.dot(Rf)) /
      std::pow(rf, 5);

  // Long range: Thole damping factors must have saturated to 1, so the
  // screened (tiny alpha) energy must match the undamped point-dipole
  // formula closely.
  BOOST_CHECK_CLOSE(e_far, e_undamped_far, 1e-3);
}

// The erf- and erfc-screened B-functions are complements by construction:
// together they must reconstruct the bare Coulomb derivatives. Asserted
// as an identity rather than against stored numbers, so it cannot be
// satisfied by a reference that mirrors a mistake in either branch.
BOOST_AUTO_TEST_CASE(erf_and_erfc_bfunctions_are_complementary) {
  EwaldRealSpaceInteractor interactor(0.3, 0.39);
  double worst = 0.0;
  for (double r : {0.5, 1.0, 2.5, 7.0, 15.0}) {
    const auto bc = interactor.ComputeB(r);
    const auto be = interactor.ComputeErfB(r);
    worst = std::max(worst, std::abs(bc.B0 + be.B0 - 1.0 / r) * r);
    worst = std::max(
        worst, std::abs(bc.B1 + be.B1 - 1.0 / std::pow(r, 3)) * std::pow(r, 3));
    worst = std::max(worst, std::abs(bc.B2 + be.B2 - 3.0 / std::pow(r, 5)) *
                                std::pow(r, 5) / 3.0);
  }
  BOOST_CHECK_SMALL(worst, 1e-12);
}

// The same complementarity, one level up: the erfc-screened induced field
// plus the erf-screened correction must equal the BARE dipole field.
//
// This is the check that pins the correction's sign convention without
// re-deriving it. A sign error in ApplyErfInducedFieldCorrection cannot
// pass this, whereas comparing against a hand-written expected vector
// could -- if that vector were derived the same (wrong) way.
//
// The separation is chosen so Thole damping is inactive (au3 >= 40),
// because the erf correction is deliberately undamped: the
// reciprocal-space contribution it removes has no damping of its own.
BOOST_AUTO_TEST_CASE(erf_correction_completes_the_bare_induced_field) {
  const double alpha = 0.3;
  const double thole = 0.39;
  EwaldRealSpaceInteractor interactor(alpha, thole);

  auto make_site = [](const Eigen::Vector3d& pos) {
    PolarSite site(0, "C", pos);
    site.setpolarization(2.0 * Eigen::Matrix3d::Identity());
    return site;
  };

  PolarSite source = make_site(Eigen::Vector3d::Zero());
  const Eigen::Vector3d mu(3e-3, -1e-3, 2e-3);
  source.setInduced_Dipole(mu);
  const Eigen::Vector3d target_pos(6.0, 2.0, -1.0);

  PolarSite screened = make_site(target_pos);
  screened.Reset();
  interactor.ApplyInducedField<Estatic::V>(source, screened);

  PolarSite correction = make_site(target_pos);
  correction.Reset();
  interactor.ApplyErfInducedFieldCorrection<Estatic::V>(source, correction);

  const double r = target_pos.norm();
  const double au3 = thole * std::pow(r, 3) * source.getSqrtInvEigenDamp() *
                     screened.getSqrtInvEigenDamp();
  BOOST_REQUIRE_GT(au3, 40.0);  // Thole must be inactive for this identity

  const Eigen::Vector3d bare =
      mu.dot(target_pos) * (3.0 / std::pow(r, 5)) * target_pos -
      (1.0 / std::pow(r, 3)) * mu;

  // The correction is stored negated (it is subtracted), so recombining
  // means screened MINUS correction.
  const Eigen::Vector3d recombined = screened.V() - correction.V();
  BOOST_CHECK_SMALL((recombined - bare).norm() / bare.norm(), 1e-12);
}

// Passing a periodic-image shift must be identical to physically moving
// the source, for the correction as much as for the screened field.
BOOST_AUTO_TEST_CASE(erf_correction_shift_matches_moving_the_source) {
  EwaldRealSpaceInteractor interactor(0.3, 0.39);
  const Eigen::Vector3d mu(3e-3, -1e-3, 2e-3);
  const Eigen::Vector3d target_pos(6.0, 2.0, -1.0);

  auto make_site = [](const Eigen::Vector3d& pos) {
    PolarSite site(0, "C", pos);
    site.setpolarization(2.0 * Eigen::Matrix3d::Identity());
    return site;
  };

  PolarSite at_origin = make_site(Eigen::Vector3d::Zero());
  at_origin.setInduced_Dipole(mu);
  PolarSite moved = make_site(Eigen::Vector3d(-4.0, 1.5, 3.0));
  moved.setInduced_Dipole(mu);
  const Eigen::Vector3d shift(4.0, -1.5, -3.0);

  PolarSite a = make_site(target_pos);
  a.Reset();
  interactor.ApplyErfInducedFieldCorrection<Estatic::V>(at_origin, a);
  PolarSite b = make_site(target_pos);
  b.Reset();
  interactor.ApplyErfInducedFieldCorrection<Estatic::V>(moved, b, shift);

  BOOST_CHECK_SMALL((a.V() - b.V()).norm() / a.V().norm(), 1e-14);
}

// The energy counterpart of erf_correction_completes_the_bare_induced_field:
// the erfc- and erf-screened permanent-multipole energies must sum to the
// bare Coulomb energy. The reference is derived here from the multipole
// expansion directly, NOT from the interactor, so it cannot agree with a
// mistake in either screened branch.
BOOST_AUTO_TEST_CASE(erfc_and_erf_static_energies_sum_to_bare) {
  EwaldRealSpaceInteractor interactor(0.3, 0.39);

  auto make = [](const Eigen::Vector3d& pos, double q,
                 const Eigen::Vector3d& mu) {
    PolarSite site(0, "C", pos);
    site.setpolarization(2.0 * Eigen::Matrix3d::Identity());
    site.setCharge(q);
    site.setStaticDipole(mu);
    return site;
  };

  PolarSite a =
      make(Eigen::Vector3d::Zero(), 0.35, Eigen::Vector3d(2e-2, -1e-2, 5e-3));
  PolarSite b = make(Eigen::Vector3d(4.0, 1.5, -2.0), -0.22,
                     Eigen::Vector3d(-8e-3, 4e-3, 1e-2));

  const double e_erfc = interactor.CalcStaticEnergy<PolarSite, PolarSite>(a, b);
  const double e_erf =
      interactor.CalcErfStaticEnergy<PolarSite, PolarSite>(a, b);

  const Eigen::Vector3d r = b.getPos() - a.getPos();
  const double R = r.norm();
  const double R3 = R * R * R;
  const double R5 = R3 * R * R;
  const Eigen::Vector3d m1 = a.getStaticDipole();
  const Eigen::Vector3d m2 = b.getStaticDipole();
  const double bare =
      a.getCharge() * b.getCharge() / R +
      (b.getCharge() * m1.dot(r) - a.getCharge() * m2.dot(r)) / R3 +
      (m1.dot(m2) / R3 - 3.0 * m1.dot(r) * m2.dot(r) / R5);

  BOOST_CHECK_SMALL(std::abs(e_erfc + e_erf - bare) / std::abs(bare), 1e-12);
}

// Coincident sites are the common case for the energy correction, not an
// edge case: every foreground target sits inside a copy being removed.
// The erf-screened energy must take its analytic r -> 0 limit there
// rather than dividing by zero.
BOOST_AUTO_TEST_CASE(erf_static_energy_is_finite_for_coincident_sites) {
  EwaldRealSpaceInteractor interactor(0.3, 0.39);
  PolarSite site(0, "C", Eigen::Vector3d::Zero());
  site.setpolarization(2.0 * Eigen::Matrix3d::Identity());
  site.setCharge(0.35);
  site.setStaticDipole(Eigen::Vector3d(2e-2, -1e-2, 5e-3));

  const double e =
      interactor.CalcErfStaticEnergy<PolarSite, PolarSite>(site, site);
  BOOST_CHECK(std::isfinite(e));

  // The charge-charge part is q^2 * 2*alpha/sqrt(pi); with a dipole also
  // present the total must still be dominated by it and positive here.
  BOOST_CHECK_GT(e, 0.0);
}

// A periodic-image shift must be equivalent to moving the source, for the
// energy as much as for the field.
BOOST_AUTO_TEST_CASE(static_energy_shift_matches_moving_the_source) {
  EwaldRealSpaceInteractor interactor(0.3, 0.39);
  auto make = [](const Eigen::Vector3d& pos) {
    PolarSite site(0, "C", pos);
    site.setpolarization(2.0 * Eigen::Matrix3d::Identity());
    site.setCharge(0.35);
    site.setStaticDipole(Eigen::Vector3d(2e-2, -1e-2, 5e-3));
    return site;
  };
  PolarSite target = make(Eigen::Vector3d(4.0, 1.5, -2.0));
  PolarSite at_origin = make(Eigen::Vector3d::Zero());
  PolarSite moved = make(Eigen::Vector3d(-7.0, 3.0, 5.0));
  const Eigen::Vector3d shift(7.0, -3.0, -5.0);

  const double direct =
      interactor.CalcStaticEnergy<PolarSite, PolarSite>(at_origin, target);
  const double shifted =
      interactor.CalcStaticEnergy<PolarSite, PolarSite>(moved, target, shift);
  BOOST_CHECK_SMALL(std::abs(shifted - direct) / std::abs(direct), 1e-14);
}

// The four erf-screened routines special-case coincident sites
// (r < kCoincidenceTol = 1e-2) with an analytic self-term. Those are
// LIMITS of the general branch, so the two must agree where they meet --
// a continuous function of r cannot jump.
//
// Not a formality: all four had the dipole term the wrong way round, a
// jump of exactly 2*self_coeff*mu. Only ApplyErfInducedFieldCorrection's
// was observable at rank 0 (the others contract against a STATIC dipole,
// zero when the .mps carry only charges), and it surfaced as an alpha^2.9
// drift in the polar energy of an 18-segment job.
//
// Tolerance is set from the actual truncation at kCoincidenceTol, not
// guessed: the general branch approaches its limit as 1.45*(alpha*r)^2
// (measured 7.12e-5, 1.78e-5, 4.40e-6 at r = 2e-2, 1e-2, 5e-3 with
// alpha = 0.35). That is 1.8e-5 here, so 1e-5 is unusable -- an earlier
// version used it and failed on correct code -- and smaller r does not
// help, since 1/r^3 - B1_erfc has lost four digits by r = 1e-4. 1e-3 is
// fifty times the truncation and still 2000x tighter than the factor of
// two a flipped sign gives.
BOOST_AUTO_TEST_CASE(coincidence_branches_match_their_own_limit) {
  const double alpha = 0.35;
  EwaldRealSpaceInteractor interactor(alpha, 0.39);
  const double self_coeff =
      (4.0 / 3.0) * alpha * alpha * alpha / std::sqrt(votca::tools::conv::Pi);

  const double r_edge = 1e-2;  // == kCoincidenceTol
  const Eigen::Vector3d mu_src(3.1e-2, -1.7e-2, 4.4e-2);
  const Eigen::Vector3d mu_tgt(-2.3e-2, 5.2e-2, 1.1e-2);

  auto make = [](const Eigen::Vector3d& pos) {
    PolarSite site(0, "C", pos);
    site.setpolarization(2.0 * Eigen::Matrix3d::Identity());
    return site;
  };

  // First, the limit itself, stated independently of any branch: the
  // general branch's dipole field approaches -self_coeff*mu, NOT +. The
  // routine applies V -= src.field, so the delivered V approaches
  // +self_coeff*mu. Getting this backwards is the whole bug.
  {
    PolarSite src = make(Eigen::Vector3d::Zero());
    src.setInduced_Dipole(mu_src);
    PolarSite tgt = make(r_edge * Eigen::Vector3d::UnitZ());
    tgt.Reset();
    interactor.ApplyErfInducedFieldCorrection<Estatic::V>(src, tgt);
    const Eigen::Vector3d expected = self_coeff * mu_src;
    std::cout << "near-coincident V = " << tgt.V().transpose()
              << "   expected = " << expected.transpose() << std::endl;
    BOOST_CHECK_SMALL((tgt.V() - expected).norm() / expected.norm(), 1e-3);
  }

  // (a) induced-dipole field correction: r = 0 against r = kCoincidenceTol
  {
    PolarSite src = make(Eigen::Vector3d::Zero());
    src.setInduced_Dipole(mu_src);

    PolarSite tgt_at = make(Eigen::Vector3d::Zero());
    tgt_at.Reset();
    interactor.ApplyErfInducedFieldCorrection<Estatic::V>(src, tgt_at);

    PolarSite tgt_near = make(r_edge * Eigen::Vector3d::UnitZ());
    tgt_near.Reset();
    interactor.ApplyErfInducedFieldCorrection<Estatic::V>(src, tgt_near);

    BOOST_CHECK_SMALL((tgt_at.V() - tgt_near.V()).norm() / tgt_near.V().norm(),
                      1e-3);
  }

  // (b) static-dipole field correction. No charge on the source: its own
  //     limit is a separate (and already correct) one, and including it
  //     would let a large right answer mask a small wrong one.
  {
    PolarSite src = make(Eigen::Vector3d::Zero());
    src.setStaticDipole(mu_src);

    PolarSite tgt_at = make(Eigen::Vector3d::Zero());
    tgt_at.Reset();
    interactor.ApplyErfStaticFieldCorrection<PolarSite, Estatic::V>(src,
                                                                    tgt_at);

    PolarSite tgt_near = make(r_edge * Eigen::Vector3d::UnitZ());
    tgt_near.Reset();
    interactor.ApplyErfStaticFieldCorrection<PolarSite, Estatic::V>(src,
                                                                    tgt_near);

    BOOST_CHECK_SMALL((tgt_at.V() - tgt_near.V()).norm() / tgt_near.V().norm(),
                      1e-3);
  }

  // (c) static-static erf energy, dipole-dipole part only (no charges,
  //     same reasoning as (b)).
  {
    PolarSite src = make(Eigen::Vector3d::Zero());
    src.setStaticDipole(mu_src);
    PolarSite tgt_at = make(Eigen::Vector3d::Zero());
    tgt_at.setStaticDipole(mu_tgt);
    PolarSite tgt_near = make(r_edge * Eigen::Vector3d::UnitZ());
    tgt_near.setStaticDipole(mu_tgt);

    const double e_at =
        interactor.CalcErfStaticEnergy<PolarSite, PolarSite>(src, tgt_at);
    const double e_near =
        interactor.CalcErfStaticEnergy<PolarSite, PolarSite>(src, tgt_near);
    std::cout << "CalcErfStaticEnergy        self = " << e_at
              << "   near = " << e_near << std::endl;
    BOOST_CHECK_SMALL(std::abs(e_at - e_near) / std::abs(e_near), 1e-3);
  }

  // (d) induced-source erf energy: source induced dipole, target static.
  {
    PolarSite src = make(Eigen::Vector3d::Zero());
    src.setInduced_Dipole(mu_src);
    PolarSite tgt_at = make(Eigen::Vector3d::Zero());
    tgt_at.setStaticDipole(mu_tgt);
    PolarSite tgt_near = make(r_edge * Eigen::Vector3d::UnitZ());
    tgt_near.setStaticDipole(mu_tgt);

    const double e_at = interactor.CalcErfInducedSourceEnergy(src, tgt_at);
    const double e_near = interactor.CalcErfInducedSourceEnergy(src, tgt_near);
    std::cout << "CalcErfInducedSourceEnergy self = " << e_at
              << "   near = " << e_near << std::endl;
    BOOST_CHECK_SMALL(std::abs(e_at - e_near) / std::abs(e_near), 1e-3);
  }
}

// What the induced real-space half must add up to with the reciprocal
// side: the Thole-damped BARE interaction, with no alpha in it.
//
// ApplyErfInducedFieldCorrection computes the undamped erf field the
// reciprocal sum contributes, stored negated since its job is to remove
// it. So applying the two routines to separate targets and SUBTRACTING
// reconstructs the pair's actual contribution:
//
//     ApplyInducedField - ApplyErfInducedFieldCorrection
//        = (l*B_bare - B_erf) - (-B_erf) = l*B_bare
//
// The reference is the textbook damped point-dipole tensor, not a
// B-function, so it cannot agree with a mistake in either branch. Under
// the old l*B_erfc convention the same combination is alpha-dependent,
// so the three alphas discriminate on their own.
BOOST_AUTO_TEST_CASE(damped_real_and_erf_halves_reassemble_to_damped_bare) {
  const double thole = 0.39;
  const Eigen::Vector3d r_vec(1.1, -1.4, 0.7);  // |r| ~ 1.93 bohr
  const Eigen::Vector3d mu(3.1e-2, -1.7e-2, 4.4e-2);

  auto make = [](const Eigen::Vector3d& pos) {
    PolarSite site(0, "C", pos);
    site.setpolarization(2.0 * Eigen::Matrix3d::Identity());
    return site;
  };

  const double r = r_vec.norm();
  Eigen::Vector3d reference;
  {
    // Damping factors, from the model directly rather than via the
    // interactor, so this reference stands on its own.
    PolarSite a = make(Eigen::Vector3d::Zero());
    PolarSite b = make(r_vec);
    const double au3 = thole * std::pow(r, 3) * a.getSqrtInvEigenDamp() *
                       b.getSqrtInvEigenDamp();
    const double expUa = std::exp(-au3);
    const double l3 = 1.0 - expUa;
    const double l5 = 1.0 - (1.0 + au3) * expUa;
    std::cout << "reassembly test: r = " << r << " bohr, l3 = " << l3
              << ", l5 = " << l5 << std::endl;
    // Guard against a vacuous pass: with l3 = l5 = 1 the old and new
    // conventions coincide and this test would prove nothing.
    BOOST_REQUIRE_LT(l3, 0.99);
    BOOST_REQUIRE_LT(l5, 0.99);

    reference = mu.dot(r_vec) * l5 * (3.0 / std::pow(r, 5)) * r_vec -
                l3 * (1.0 / std::pow(r, 3)) * mu;
  }

  std::vector<Eigen::Vector3d> assembled;
  for (double alpha : {0.20, 0.35, 0.60}) {
    EwaldRealSpaceInteractor interactor(alpha, thole);

    PolarSite source = make(Eigen::Vector3d::Zero());
    source.setInduced_Dipole(mu);

    PolarSite target_real = make(r_vec);
    target_real.Reset();
    interactor.ApplyInducedField<Estatic::V>(source, target_real);

    PolarSite target_erf = make(r_vec);
    target_erf.Reset();
    interactor.ApplyErfInducedFieldCorrection<Estatic::V>(source, target_erf);

    const Eigen::Vector3d total = target_real.V() - target_erf.V();
    std::cout << "  alpha = " << alpha
              << "  reassembled = " << total.transpose() << std::endl;
    assembled.push_back(total);

    BOOST_CHECK_SMALL((total - reference).norm() / reference.norm(), 1e-12);
  }

  std::cout << "  reference   = " << reference.transpose() << std::endl;

  // Stated separately, so that a failure of alpha-independence is
  // distinguishable from a failure to match the analytic form.
  for (std::size_t i = 1; i < assembled.size(); ++i) {
    BOOST_CHECK_SMALL(
        (assembled[i] - assembled[0]).norm() / assembled[0].norm(), 1e-12);
  }
}

BOOST_AUTO_TEST_SUITE_END()
