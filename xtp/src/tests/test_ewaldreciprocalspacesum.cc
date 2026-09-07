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

#define BOOST_TEST_MODULE ewaldreciprocalspacesum_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/ewaldreciprocalspacesum.h"
#include "votca/xtp/ewaldrealspacesum.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(ewaldreciprocalspacesum_test)

PolarSegment MakeChargeSegment(Index id, const Eigen::Vector3d& pos,
                               double charge) {
  PolarSegment seg("seg", id);
  PolarSite site(id, "H", pos);
  site.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles = Vector9d::Zero();
  mpoles(0) = charge;
  site.setMultipole(mpoles, 0);
  seg.push_back(site);
  return seg;
}

// The real test of a real-space/reciprocal-space split: the *combined*
// field (real+reciprocal together, at a fixed pair of well-converged
// cutoffs) must not depend on the arbitrary Ewald splitting parameter
// alpha, since alpha is purely a computational convenience that decides
// how much of the true 1/r interaction is pushed into real space versus
// reciprocal space. This does not require the shape/k=0 correction (out
// of scope for this class) to hold, since that term is itself alpha-
// independent -- it is a genuinely strong, sign-sensitive check on its
// own, and one that could not be performed on the real-space sum in
// isolation (see test_ewaldrealspacesum.cc), since the real-space term
// alone has no alpha-independent target to compare against.
BOOST_AUTO_TEST_CASE(alpha_independence_of_combined_field) {
  Eigen::Matrix3d box = 30.0 * Eigen::Matrix3d::Identity();

  auto compute_combined_field = [&](double alpha) {
    EwaldRegistry registry;
    registry.Register(1, EwaldChargeState::Neutral,
                      MakeChargeSegment(1, Eigen::Vector3d::Zero(), 0.6));

    PolarSite target(2, "H", Eigen::Vector3d(2.5, 1.0, -0.7));

    // k_max must be generous enough for the larger alpha value tested
    // below (larger alpha means a slower-decaying Gaussian weight in
    // k-space, i.e. more k-vectors needed for convergence).
    EwaldReciprocalSpaceSum recip(box, registry, alpha, /*k_max=*/25.0);
    recip.AddFieldAt<Estatic::V>(target, EwaldChargeState::Neutral);

    EwaldRealSpaceSum real(box, registry, alpha, 0.39, /*r_min=*/8.0,
                           /*field_tol=*/1e-12);
    real.AddFieldAt<Estatic::V>(2, target, EwaldChargeState::Neutral);

    return target.V();
  };

  Eigen::Vector3d field_a = compute_combined_field(0.3);
  Eigen::Vector3d field_b = compute_combined_field(0.5);

  BOOST_CHECK(field_a.isApprox(field_b, 1e-4));
  if (!field_a.isApprox(field_b, 1e-4)) {
    std::cout << "alpha=0.3: " << field_a.transpose() << std::endl;
    std::cout << "alpha=0.5: " << field_b.transpose() << std::endl;
  }
  // Sanity: not vacuously passing because both are zero.
  BOOST_CHECK(field_a.norm() > 1e-3);
}

// A charge on the x-axis, target also on the x-axis: by symmetry, the
// reciprocal-space field at the target must lie entirely along x (y and z
// components must vanish), and -- same sign convention as a repulsive
// Coulomb field from a positive source -- must point away from the
// source charge along that axis, matching EwaldRealSpaceInteractor's own
// (already-validated) sign convention.
BOOST_AUTO_TEST_CASE(onaxis_symmetry_and_sign) {
  Eigen::Matrix3d box = 10.0 * Eigen::Matrix3d::Identity();
  EwaldRegistry registry;
  registry.Register(1, EwaldChargeState::Neutral,
                    MakeChargeSegment(1, Eigen::Vector3d::Zero(), 1.0));

  PolarSite target(2, "H", Eigen::Vector3d(3.0, 0.0, 0.0));

  EwaldReciprocalSpaceSum recip(box, registry, 0.35, /*k_max=*/20.0);
  recip.AddFieldAt<Estatic::V>(target, EwaldChargeState::Neutral);

  BOOST_CHECK_SMALL(target.V().y(), 1e-10);
  BOOST_CHECK_SMALL(target.V().z(), 1e-10);
  BOOST_CHECK(target.V().x() > 0.0);
}

// EwaldReciprocalSpaceSum, unlike EwaldRealSpaceSum, never excludes
// anything from its own global structure factor -- not even the
// target's own registered segment (see this class's own class
// documentation for why: legacy PolarBackground's own reciprocal-space
// code, both for the permanent case (SP_SFactorCalc/FP_KFieldCalc) and
// the induced case (SU_SFactorCalc/FU_KFieldCalc), works the same way,
// and the Ewald-split identity erfc(a*r)/r + [reciprocal-space
// contribution] = 1/r requires this reciprocal-space contribution to be
// present, unconditionally, for the two halves to combine correctly
// wherever EwaldRealSpaceSum's own real-space sum legitimately excludes
// a pair). An earlier version of this test checked the OPPOSITE property
// (intramolecular contribution excluded) -- correct for the code as it
// existed then, but wrong now.
//
// What this test actually checks: a target DOES feel a nonzero
// reciprocal-space contribution from another site of its own registered
// segment -- but not a nonzero contribution from its own literal self.
// The self-contribution is zero by construction, not by exclusion: a
// site's own term in the global structure factor sum is q*exp(-i*k.r)
// (purely real when evaluated back at that same site's own position r,
// since the phase exactly cancels), and Im[V] of a purely real quantity
// is 0 -- so a point charge/dipole never feels a net field from itself
// via this formula, regardless of whether anything is excluded.
BOOST_AUTO_TEST_CASE(intramolecular_contribution_included) {
  Eigen::Matrix3d box = 10.0 * Eigen::Matrix3d::Identity();
  EwaldRegistry registry;
  PolarSegment seg("seg", 1);
  PolarSite site_a(1, "H", Eigen::Vector3d::Zero());
  site_a.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles_a = Vector9d::Zero();
  mpoles_a(0) = 1.0;
  site_a.setMultipole(mpoles_a, 0);
  seg.push_back(site_a);

  PolarSite site_b(2, "H", Eigen::Vector3d(2.0, 0.0, 0.0));
  site_b.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles_b = Vector9d::Zero();
  mpoles_b(0) = -1.0;
  site_b.setMultipole(mpoles_b, 0);
  seg.push_back(site_b);

  registry.Register(1, EwaldChargeState::Neutral, seg);

  PolarSite target = site_b;
  target.Reset();

  EwaldReciprocalSpaceSum recip(box, registry, 0.35, /*k_max=*/20.0);
  recip.AddFieldAt<Estatic::V>(target, EwaldChargeState::Neutral);

  // Not zero: site_b now genuinely feels site_a's own contribution
  // (opposite charges 2.0 bohr apart -- a real, sizeable field, not a
  // rounding-level residual).
  BOOST_CHECK(target.V().norm() > 1e-3);

  // Cross-check against an INDEPENDENT registry containing only site_a
  // (no site_b at all) -- the field a lone external charge would
  // produce at the same target position, entirely unrelated to any
  // exclusion/inclusion question about site_b's own segment membership.
  // If intramolecular contributions are genuinely included (not merely
  // "not exactly excluded" through some other accident), this must match
  // the two-site registry's own result exactly, since site_b's own
  // self-contribution is exactly zero (see class-level comment above).
  EwaldRegistry registry_external_only;
  PolarSegment seg_a_only("seg", 1);
  seg_a_only.push_back(site_a);
  registry_external_only.Register(1, EwaldChargeState::Neutral,
                                  seg_a_only);
  PolarSite target_external = site_b;
  target_external.Reset();
  EwaldReciprocalSpaceSum recip_external(box, registry_external_only, 0.35,
                                         /*k_max=*/20.0);
  recip_external.AddFieldAt<Estatic::V>(target_external,
                                        EwaldChargeState::Neutral);

  BOOST_CHECK(target.V().isApprox(target_external.V(), 1e-12));
}

BOOST_AUTO_TEST_SUITE_END()
