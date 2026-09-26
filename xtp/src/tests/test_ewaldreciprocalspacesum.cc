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
#include <memory>
#include <vector>

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/ewaldrealspaceinteractor.h"
#include "votca/xtp/ewaldrealspacesum.h"
#include "votca/xtp/ewaldreciprocalspacesum.h"
#include "votca/xtp/ewaldshapecorrection.h"

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

    // k_max must scale with alpha, because a larger alpha means a
    // slower-decaying Gaussian weight in k-space. Tying it to alpha as
    // 12*alpha rather than picking one flat literal generous enough for
    // the largest alpha keeps the truncation error FIXED across the
    // scan -- exp(-(12*alpha)^2/(4*alpha^2)) = exp(-36) = 2.3e-16,
    // against a 1e-4 comparison below -- which is what this test needs:
    // an alpha-dependent truncation error would imitate the failure it
    // is looking for. It is also far cheaper, cubically so: 98,384
    // k-vectors at alpha = 0.5 in this 30 bohr box, against 7,123,810
    // at the flat k_max = 25.0 this used to pass.
    EwaldReciprocalSpaceSum recip(box, registry, alpha,
                                  /*k_max=*/12.0 * alpha);
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
  registry_external_only.Register(1, EwaldChargeState::Neutral, seg_a_only);
  PolarSite target_external = site_b;
  target_external.Reset();
  EwaldReciprocalSpaceSum recip_external(box, registry_external_only, 0.35,
                                         /*k_max=*/20.0);
  recip_external.AddFieldAt<Estatic::V>(target_external,
                                        EwaldChargeState::Neutral);

  BOOST_CHECK(target.V().isApprox(target_external.V(), 1e-12));
}

// THE INVARIANT THAT VALIDATES THE WHOLE SPLITTING.
//
// alpha decides only how the Coulomb interaction is divided between the
// real- and reciprocal-space sums. It has no physical meaning, so the
// TOTAL must not depend on it, even though each half changes strongly.
//
// This is far more sensitive than checking either half alone, and it is
// what caught a real error while this was being written: an erf
// correction had been subtracted from the energy by analogy with the
// FIELD, where it is genuinely needed because the field's reciprocal
// part sums the total structure factor. The energy's reciprocal part is
// a CROSS term that never contains the foreground, so the correction
// removed a foreground-foreground contribution that was never added.
// Being erf-weighted, the spurious term is alpha-dependent, so it broke
// this invariant immediately while leaving each half individually
// plausible.
BOOST_AUTO_TEST_CASE(static_energy_is_independent_of_the_splitting) {
  const double L = 20.0;
  const double thole = 0.39;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;
  const Index n = 3;
  const double d = L / double(n);
  Index id = 0;
  for (Index a = 0; a < n; ++a) {
    for (Index b = 0; b < n; ++b) {
      for (Index c = 0; c < n; ++c) {
        const Eigen::Vector3d centre(double(a) * d, double(b) * d,
                                     double(c) * d);
        PolarSegment seg("seg", id);
        const double t = 0.63;
        const Eigen::Vector3d offsets[3] = {
            {0.0, 0.0, 0.0}, {t, t, t}, {-t, -t, t}};
        const double charges[3] = {-0.4, 0.2, 0.2};
        for (Index j = 0; j < 3; ++j) {
          PolarSite site(j, (j == 0) ? "C" : "H", centre + offsets[j]);
          site.setpolarization(3.0 * Eigen::Matrix3d::Identity());
          site.setCharge(charges[j]);
          site.setStaticDipole(
              Eigen::Vector3d(5e-3 * double(j + 1), -2e-3, 1e-3 * double(j)));
          seg.push_back(site);
        }
        registry.Register(id, EwaldChargeState::Neutral, seg);
        ++id;
      }
    }
  }

  const Index fg_id = 0;
  const PolarSegment& fg_seg = registry.Get(fg_id, EwaldChargeState::Neutral);
  Eigen::Vector3d fg_pos = Eigen::Vector3d::Zero();
  Index n_sites = 0;
  for (const PolarSite& site : fg_seg) {
    fg_pos += site.getPos();
    ++n_sites;
  }
  fg_pos /= double(n_sites);
  const std::vector<std::pair<Index, Eigen::Vector3d>> foreground{
      {fg_id, fg_pos}};

  double reference = 0.0;
  bool have_reference = false;
  for (double alpha : {0.20, 0.25, 0.30, 0.35}) {
    // Both cutoffs are scaled with alpha so each sum stays equally well
    // converged; otherwise this would measure truncation, not the
    // splitting.
    const double k_max = 14.0 * alpha;
    const double r_min = 6.0 / alpha;

    EwaldRealSpaceSum real_sum(box, registry, alpha, thole, r_min, 1e-14, 0.945,
                               25, 6.0, foreground);
    EwaldReciprocalSpaceSum recip_sum(box, registry, alpha, k_max);

    double e_real = 0.0;
    std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> fg_sites;
    for (const PolarSite& site : fg_seg) {
      PolarSite probe = site;
      probe.Reset();
      real_sum.AddFieldAt<Estatic::V>(fg_id, probe, EwaldChargeState::Neutral);
      e_real += real_sum.CalcStaticEnergyAt(probe, EwaldChargeState::Neutral);
      fg_sites.push_back({&site, site.getPos()});
    }
    // NOTHING is held out of S_bg (this revision). The foreground
    // segment's COINCIDENT copy is what real space suppresses, and it is
    // removed from the reciprocal side by subtracting its erf energy --
    // its other periodic images are ordinary background molecules and
    // belong in both sums. Holding the segment out of S_bg instead
    // removed those images from the reciprocal side while real space
    // kept them, which is exactly the alpha-dependence this case
    // measures. See EwaldRegion's step (5).
    const std::vector<const PolarSite*> no_exclusions;
    const double e_recip = recip_sum.CalcStaticEnergyBetween(
        fg_sites, no_exclusions, EwaldChargeState::Neutral);

    EwaldRealSpaceInteractor inter(alpha, thole);
    double e_erf = 0.0;
    for (const PolarSite& source : fg_seg) {
      for (const auto& entry : fg_sites) {
        e_erf += inter.CalcErfStaticEnergy<PolarSite, PolarSite>(
            source, *entry.first, Eigen::Vector3d::Zero());
      }
    }
    const double total = e_real + e_recip - e_erf;

    if (!have_reference) {
      reference = total;
      have_reference = true;
      // Guard against the invariant holding trivially because everything
      // is zero.
      BOOST_REQUIRE_GT(std::abs(total), 1e-12);
    } else {
      BOOST_CHECK_SMALL(std::abs(total - reference) / std::abs(reference),
                        1e-8);
    }
  }
}

// A 3x3x3 lattice of three-site segments carrying BOTH permanent
// charges and induced dipoles, for the induced-source tests below. The
// induced dipoles are deliberately unequal and non-symmetric so that no
// component of the term can cancel by accident; the permanent charges
// sum to zero per segment, keeping the cell neutral so the omitted k=0
// term is exactly zero.
EwaldRegistry BuildInducedTestRegistry(double L) {
  EwaldRegistry registry;
  const Index n = 3;
  const double d = L / double(n);
  Index id = 0;
  for (Index a = 0; a < n; ++a) {
    for (Index b = 0; b < n; ++b) {
      for (Index c = 0; c < n; ++c) {
        const Eigen::Vector3d centre(double(a) * d, double(b) * d,
                                     double(c) * d);
        PolarSegment seg("seg", id);
        const double t = 0.63;
        const Eigen::Vector3d offsets[3] = {
            {0.0, 0.0, 0.0}, {t, t, t}, {-t, -t, t}};
        const double charges[3] = {-0.4, 0.2, 0.2};
        for (Index j = 0; j < 3; ++j) {
          PolarSite site(j, (j == 0) ? "C" : "H", centre + offsets[j]);
          site.setpolarization(3.0 * Eigen::Matrix3d::Identity());
          site.setCharge(charges[j]);
          // Induced dipoles vary with both the site index and the cell
          // index, so the background's total induced moment does not
          // vanish by symmetry.
          site.setInduced_Dipole(Eigen::Vector3d(
              1e-2 * double(j + 1) + 1e-3 * double(a), -7e-3 + 2e-3 * double(b),
              4e-3 * double(j) - 1e-3 * double(c)));
          seg.push_back(site);
        }
        registry.Register(id, EwaldChargeState::Neutral, seg);
        ++id;
      }
    }
  }
  return registry;
}

// The same invariant for the [foreground permanent] x [background
// induced] term, which is split across real and reciprocal space the
// same way -- so the same alpha-independence must hold, and for the same
// reason.
//
// With ONE caveat that shapes this test: the real-space half of that
// term is Thole-damped and the reciprocal half is not (legacy has the
// identical asymmetry -- FU12_ERFC_At_By damps, FU12_ERF_At_By does
// not). Damping therefore leaks alpha-dependence into the total,
// because changing alpha moves weight between a damped sum and an
// undamped one in exactly the region where damping is active.
//
// So the invariant is tested with damping switched OFF, which is what a
// very large thole_a does: ComputeThole leaves l3 = l5 = 1 whenever
// au3 >= 40. That makes the decomposition exact, and the test then
// genuinely checks the real/reciprocal split rather than the damping
// model. A second case below measures how large the damped
// alpha-dependence actually is, so the approximation is quantified
// rather than assumed small.
BOOST_AUTO_TEST_CASE(induced_source_energy_is_independent_of_the_splitting) {
  const double L = 20.0;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  // Large enough that au3 >= 40 for every pair that matters, i.e. no
  // damping anywhere -- see this test's own comment above.
  const double thole_undamped = 1e6;

  EwaldRegistry registry = BuildInducedTestRegistry(L);

  const Index fg_id = 0;
  const PolarSegment& fg_seg = registry.Get(fg_id, EwaldChargeState::Neutral);
  Eigen::Vector3d fg_pos = Eigen::Vector3d::Zero();
  Index n_sites = 0;
  for (const PolarSite& site : fg_seg) {
    fg_pos += site.getPos();
    ++n_sites;
  }
  fg_pos /= double(n_sites);
  const std::vector<std::pair<Index, Eigen::Vector3d>> foreground{
      {fg_id, fg_pos}};

  double reference = 0.0;
  bool have_reference = false;
  for (double alpha : {0.20, 0.25, 0.30, 0.35}) {
    const double k_max = 14.0 * alpha;
    const double r_min = 6.0 / alpha;

    EwaldRealSpaceSum real_sum(box, registry, alpha, thole_undamped, r_min,
                               1e-14, 0.945, 25, 6.0, foreground);
    EwaldReciprocalSpaceSum recip_sum(box, registry, alpha, k_max);

    double e_real = 0.0;
    std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> fg_sites;
    for (const PolarSite& site : fg_seg) {
      PolarSite probe = site;
      probe.Reset();
      real_sum.AddFieldAt<Estatic::V>(fg_id, probe, EwaldChargeState::Neutral);
      e_real +=
          real_sum.CalcInducedSourceEnergyAt(probe, EwaldChargeState::Neutral);
      fg_sites.push_back({&site, site.getPos()});
    }
    // Nothing held out of S_bg; the coincident copy comes off via its
    // erf energy instead. Same reasoning as the static case above.
    const std::vector<const PolarSite*> no_exclusions;
    const double e_recip = recip_sum.CalcInducedSourceEnergyBetween(
        fg_sites, no_exclusions, EwaldChargeState::Neutral);

    EwaldRealSpaceInteractor inter(alpha, thole_undamped);
    double e_erf = 0.0;
    for (const PolarSite& source : fg_seg) {
      for (const auto& entry : fg_sites) {
        e_erf += inter.CalcErfInducedSourceEnergy(source, *entry.first,
                                                  Eigen::Vector3d::Zero());
      }
    }
    const double total = e_real + e_recip - e_erf;

    if (!have_reference) {
      reference = total;
      have_reference = true;
      BOOST_REQUIRE_GT(std::abs(total), 1e-12);
    } else {
      BOOST_CHECK_SMALL(std::abs(total - reference) / std::abs(reference),
                        1e-8);
    }
  }
}

// Two things about the Thole damping on this term, both measured rather
// than assumed.
//
// The real-space half of the induced-source term is damped and the
// reciprocal half is not (legacy has the same asymmetry:
// FU12_ERFC_At_By damps, FU12_ERF_At_By does not), so in principle the
// split is only approximately alpha-independent. In practice, at
// realistic damping it is not approximate at all -- damping is active
// only at sub-molecular separations, and the real-space sum already
// excludes same-segment pairs, so every pair that reaches this term is
// undamped. Measured on this fixture: thole_a = 0.39 gives the same
// energy as thole_a = 1e6 to nine digits, with an alpha drift of 7e-10
// against 1.6e-14 for the strictly undamped case.
//
// So the first case below records that the caveat costs nothing at
// realistic parameters -- and the second makes sure damping is
// nevertheless WIRED IN, by turning it up until it bites. Without that
// second case, code that silently never applied damping would pass
// every other test in this file.
BOOST_AUTO_TEST_CASE(damping_is_inactive_at_intermolecular_range) {
  const double L = 20.0;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry = BuildInducedTestRegistry(L);

  const Index fg_id = 0;
  const PolarSegment& fg_seg = registry.Get(fg_id, EwaldChargeState::Neutral);
  Eigen::Vector3d fg_pos = Eigen::Vector3d::Zero();
  Index n_sites = 0;
  for (const PolarSite& site : fg_seg) {
    fg_pos += site.getPos();
    ++n_sites;
  }
  fg_pos /= double(n_sites);
  const std::vector<std::pair<Index, Eigen::Vector3d>> foreground{
      {fg_id, fg_pos}};

  auto total_at = [&](double alpha, double thole) {
    const double k_max = 14.0 * alpha;
    const double r_min = 6.0 / alpha;
    EwaldRealSpaceSum real_sum(box, registry, alpha, thole, r_min, 1e-14, 0.945,
                               25, 6.0, foreground);
    EwaldReciprocalSpaceSum recip_sum(box, registry, alpha, k_max);

    double e_real = 0.0;
    std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> fg_sites;
    for (const PolarSite& site : fg_seg) {
      PolarSite probe = site;
      probe.Reset();
      real_sum.AddFieldAt<Estatic::V>(fg_id, probe, EwaldChargeState::Neutral);
      e_real +=
          real_sum.CalcInducedSourceEnergyAt(probe, EwaldChargeState::Neutral);
      fg_sites.push_back({&site, site.getPos()});
    }
    // Nothing held out of S_bg; the coincident copy comes off via its
    // erf energy instead. Same reasoning as the two cases above.
    const std::vector<const PolarSite*> no_exclusions;
    const double e_recip = recip_sum.CalcInducedSourceEnergyBetween(
        fg_sites, no_exclusions, EwaldChargeState::Neutral);

    EwaldRealSpaceInteractor inter(alpha, thole);
    double e_erf = 0.0;
    for (const PolarSite& source : fg_seg) {
      for (const auto& entry : fg_sites) {
        e_erf += inter.CalcErfInducedSourceEnergy(source, *entry.first,
                                                  Eigen::Vector3d::Zero());
      }
    }
    return e_real + e_recip - e_erf;
  };

  // Realistic damping reproduces the undamped answer.
  const double undamped = total_at(0.25, 1e6);
  const double realistic = total_at(0.25, 0.39);
  BOOST_REQUIRE_GT(std::abs(undamped), 1e-12);
  BOOST_CHECK_SMALL(std::abs(realistic - undamped) / std::abs(undamped), 1e-8);

  // And the split stays alpha-independent with damping switched on.
  const double lo = total_at(0.20, 0.39);
  const double hi = total_at(0.35, 0.39);
  BOOST_CHECK_SMALL(std::abs(hi - lo) / std::abs(lo), 1e-8);
}

// Turn the damping up until it bites, to prove it is actually applied.
// At thole_a = 0.01 on this fixture the term changes sign, so this
// cannot pass on code that ignores damping entirely.
BOOST_AUTO_TEST_CASE(damping_is_applied_when_it_is_active) {
  const double L = 20.0;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry = BuildInducedTestRegistry(L);

  const Index fg_id = 0;
  const PolarSegment& fg_seg = registry.Get(fg_id, EwaldChargeState::Neutral);
  Eigen::Vector3d fg_pos = Eigen::Vector3d::Zero();
  Index n_sites = 0;
  for (const PolarSite& site : fg_seg) {
    fg_pos += site.getPos();
    ++n_sites;
  }
  fg_pos /= double(n_sites);
  const std::vector<std::pair<Index, Eigen::Vector3d>> foreground{
      {fg_id, fg_pos}};

  auto real_part_at = [&](double thole) {
    const double alpha = 0.25;
    EwaldRealSpaceSum real_sum(box, registry, alpha, thole, 6.0 / alpha, 1e-14,
                               0.945, 25, 6.0, foreground);
    double e_real = 0.0;
    for (const PolarSite& site : fg_seg) {
      PolarSite probe = site;
      probe.Reset();
      real_sum.AddFieldAt<Estatic::V>(fg_id, probe, EwaldChargeState::Neutral);
      e_real +=
          real_sum.CalcInducedSourceEnergyAt(probe, EwaldChargeState::Neutral);
    }
    return e_real;
  };

  const double undamped = real_part_at(1e6);
  const double damped = real_part_at(0.01);
  BOOST_REQUIRE_GT(std::abs(undamped), 1e-12);
  // Not a marginal difference: damping must move this substantially.
  BOOST_CHECK_GT(std::abs(damped - undamped) / std::abs(undamped), 0.1);
}

EwaldRegistry BuildNeutralLattice(double L) {
  EwaldRegistry registry;
  const Index n = 3;
  const double d = L / double(n);
  Index id = 0;
  for (Index a = 0; a < n; ++a) {
    for (Index b = 0; b < n; ++b) {
      for (Index c = 0; c < n; ++c) {
        const Eigen::Vector3d centre(double(a) * d, double(b) * d,
                                     double(c) * d);
        PolarSegment seg("seg", id);
        const double t = 0.63;
        const Eigen::Vector3d offsets[3] = {
            {0.0, 0.0, 0.0}, {t, t, t}, {-t, -t, t}};
        const double charges[3] = {-0.4, 0.2, 0.2};
        for (Index j = 0; j < 3; ++j) {
          PolarSite site(j, (j == 0) ? "C" : "H", centre + offsets[j]);
          site.setpolarization(3.0 * Eigen::Matrix3d::Identity());
          site.setCharge(charges[j]);
          site.setStaticDipole(
              Eigen::Vector3d(5e-3 * double(j + 1), -2e-3, 1e-3 * double(j)));
          seg.push_back(site);
        }
        registry.Register(id, EwaldChargeState::Neutral, seg);
        ++id;
      }
    }
  }
  return registry;
}

// The invariant the rest of this suite misses: a foreground that is BOTH
// charged and multi-segment. Charge alone passes (n_fg = 1 below is
// exact), multiple segments alone passes; together they caught a version
// that held the whole foreground out of S_bg and so deleted its periodic
// images -- a lattice of vacancies rather than one carved-out cavity,
// worth 1e-3 eV of alpha-dependence on a production job.
BOOST_AUTO_TEST_CASE(charged_multisegment_foreground_splitting) {
  const double L = 20.0;
  const double thole = 0.39;
  const double V = L * L * L;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry = BuildNeutralLattice(L);

  for (std::vector<Index> ids :
       std::vector<std::vector<Index>>{{0}, {0, 1}, {0, 1, 3, 9}}) {
    // Foreground: copies of those segments, +1 added to the first.
    std::vector<PolarSegment> fg;
    bool first = true;
    for (Index id : ids) {
      const PolarSegment& bg = registry.Get(id, EwaldChargeState::Neutral);
      PolarSegment seg("seg", id);
      bool first_site = true;
      for (const PolarSite& site : bg) {
        PolarSite copy = site;
        if (first && first_site) {
          copy.setCharge(site.getCharge() + 1.0);
          copy.setStaticDipole(site.getStaticDipole());
        }
        seg.push_back(copy);
        first_site = false;
      }
      fg.push_back(seg);
      first = false;
    }

    std::vector<std::pair<Index, Eigen::Vector3d>> foreground;
    for (const PolarSegment& f : fg) {
      const PolarSegment& bg =
          registry.Get(f.getId(), EwaldChargeState::Neutral);
      Eigen::Vector3d ctr = Eigen::Vector3d::Zero();
      Index n = 0;
      for (const PolarSite& s : bg) {
        ctr += s.getPos();
        ++n;
      }
      foreground.push_back({f.getId(), ctr / double(n)});
    }

    double reference = 0.0;
    bool have_reference = false;
    // Reaches down to 0.12 deliberately: at alpha >= 0.25 in this box the
    // erfc-screened share of a segment's own-image interaction is already
    // negligible, so a scan starting there cannot see the real-space half
    // of the split at all.
    for (double alpha : {0.12, 0.20, 0.30}) {
      const double k_max = 14.0 * alpha;
      EwaldRealSpaceSum real_sum(box, registry, alpha, thole, 6.0 / alpha,
                                 1e-14, 0.945, 40, 6.0, foreground);
      EwaldReciprocalSpaceSum recip(box, registry, alpha, k_max);
      EwaldShapeCorrection shape(V, registry, EwaldShape::Cube);
      EwaldRealSpaceInteractor inter(alpha, thole);

      double total = 0.0;
      for (const PolarSegment& f : fg) {
        for (const PolarSite& site : f) {
          PolarSite probe = site;
          probe.Reset();
          real_sum.AddFieldAt<Estatic::V>(f.getId(), probe,
                                          EwaldChargeState::Neutral);
          total +=
              real_sum.CalcStaticEnergyAt(probe, EwaldChargeState::Neutral);
        }
      }
      // Nothing held out of S_bg; every foreground copy, the target's own
      // included, comes off via its erf energy. Mirrors
      // EwaldRegion::ApplyFieldTo step (5).
      std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> si;
      for (const PolarSegment& f : fg) {
        for (const PolarSite& s : f) {
          si.push_back({&s, s.getPos()});
        }
      }
      const std::vector<const PolarSite*> no_exclusions;
      total += recip.CalcStaticEnergyBetween(si, no_exclusions,
                                             EwaldChargeState::Neutral);
      total += shape.CalcStaticEnergyBetween(si, no_exclusions,
                                             EwaldChargeState::Neutral);
      for (const PolarSegment& f : fg) {
        const PolarSegment& cp =
            registry.Get(f.getId(), EwaldChargeState::Neutral);
        for (const PolarSite& src : cp) {
          for (const auto& tgt : si) {
            total -= inter.CalcErfStaticEnergy<PolarSite, PolarSite>(
                src, *tgt.first);
          }
        }
      }

      if (!have_reference) {
        reference = total;
        have_reference = true;
        BOOST_REQUIRE_GT(std::abs(total), 1e-12);
      } else {
        // 1e-5 rather than the 1e-8 the single-segment cases reach. The
        // residual shrinks ~30x per alpha step, so it is truncation, not
        // a missing term -- but it has never been explained, and the
        // loose tolerance records that rather than hiding it.
        BOOST_CHECK_SMALL(std::abs(total - reference) / std::abs(reference),
                          1e-5);
      }
    }
  }
}

// ORACLE FOR EwaldRegion::PotentialAt.
//
// The potential at a point is what a UNIT TEST CHARGE placed there would
// report as its own energy: with q = 1 and mu = 0, every energy routine
// in this machinery reduces to q*phi - mu.E = phi. That is not a
// convenience -- it is what pins phi's gauge. phi carries an arbitrary
// additive constant (the k = 0 term is omitted, i.e. the uniform
// neutralising background), and a charged region embedded in it shifts
// by q * phi_0. Defining phi through the same expressions the energy
// already uses fixes that constant to the one every validated number in
// this code was computed with, instead of leaving it to be rediscovered
// as an unexplained offset once a QM density is sitting in the hole.
//
// So this case asserts the round trip: phi and E sampled at each
// foreground site with a unit probe, recombined as sum_i [q_i phi_i -
// mu_i . E_i], must reproduce the energy assembled directly from those
// sites. Two different entry points into the same sums; agreement to
// round-off or one of them is wrong.
//
// PROBES MUST HAVE DISTINCT, STABLE ADDRESSES. EwaldRealSpaceSum keys
// its neighbour cache on the target's address and writes an entry on
// first use, so a probe built on the stack inside a loop -- reusing one
// address for a succession of different positions -- gets the first
// point's neighbour list for every later point, silently. Held in
// unique_ptrs here for exactly that reason; the same constraint is why
// PotentialAt cannot be implemented by walking a DFT grid with a probe.
BOOST_AUTO_TEST_CASE(unit_probe_potential_reproduces_the_static_energy) {
  const double L = 20.0;
  const double thole = 0.39;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();

  // The lattice of static_energy_is_independent_of_the_splitting, but
  // with dipoles at the magnitude a real rank-1 .mps carries, so the
  // dipole half of phi is actually exercised.
  EwaldRegistry registry;
  const Index n = 3;
  const double d = L / double(n);
  Index id = 0;
  for (Index a = 0; a < n; ++a) {
    for (Index b = 0; b < n; ++b) {
      for (Index c = 0; c < n; ++c) {
        const Eigen::Vector3d centre(double(a) * d, double(b) * d,
                                     double(c) * d);
        PolarSegment seg("seg", id);
        const double t = 0.63;
        const Eigen::Vector3d offsets[3] = {
            {0.0, 0.0, 0.0}, {t, t, t}, {-t, -t, t}};
        const double charges[3] = {-0.4, 0.2, 0.2};
        for (Index j = 0; j < 3; ++j) {
          PolarSite site(j, (j == 0) ? "C" : "H", centre + offsets[j]);
          site.setpolarization(3.0 * Eigen::Matrix3d::Identity());
          site.setCharge(charges[j]);
          site.setStaticDipole(Eigen::Vector3d(
              0.15 * double(j + 1) + 1e-2 * double(a), -0.10 + 1e-2 * double(b),
              0.05 * double(j) - 1e-2 * double(c)));
          seg.push_back(site);
        }
        registry.Register(id, EwaldChargeState::Neutral, seg);
        ++id;
      }
    }
  }

  const Index fg_id = 0;
  const PolarSegment& fg_seg = registry.Get(fg_id, EwaldChargeState::Neutral);
  Eigen::Vector3d fg_pos = Eigen::Vector3d::Zero();
  Index n_sites = 0;
  for (const PolarSite& site : fg_seg) {
    fg_pos += site.getPos();
    ++n_sites;
  }
  fg_pos /= double(n_sites);
  const std::vector<std::pair<Index, Eigen::Vector3d>> foreground{
      {fg_id, fg_pos}};

  // Well away from the foreground segment, where phi has no business
  // depending on the splitting either.
  const Eigen::Vector3d far_point(0.5 * d, 0.37 * d, 0.71 * d);
  double phi_far_reference = 0.0;
  bool have_far = false;

  for (double alpha : {0.20, 0.25, 0.30, 0.35}) {
    const double k_max = 14.0 * alpha;
    const double r_min = 6.0 / alpha;

    EwaldRealSpaceSum real_sum(box, registry, alpha, thole, r_min, 1e-14, 0.945,
                               25, 6.0, foreground);
    EwaldReciprocalSpaceSum recip_sum(box, registry, alpha, k_max);
    EwaldRealSpaceInteractor inter(alpha, thole);
    const std::vector<const PolarSite*> no_exclusions;

    // ---- reference: assembled from the real foreground sites, exactly
    // as EwaldRegion does it (shape omitted, as in the sibling case --
    // it is alpha-independent and common to both sides here).
    std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> fg_sites;
    std::vector<std::unique_ptr<PolarSite>> site_probes;
    site_probes.reserve(std::size_t(n_sites));
    double e_real = 0.0;
    for (const PolarSite& site : fg_seg) {
      site_probes.push_back(std::make_unique<PolarSite>(site));
      PolarSite& probe = *site_probes.back();
      probe.Reset();
      real_sum.AddFieldAt<Estatic::V>(fg_id, probe, EwaldChargeState::Neutral);
      e_real += real_sum.CalcStaticEnergyAt(probe, EwaldChargeState::Neutral);
      fg_sites.push_back({&site, site.getPos()});
    }
    const double e_recip = recip_sum.CalcStaticEnergyBetween(
        fg_sites, no_exclusions, EwaldChargeState::Neutral);
    double e_erf = 0.0;
    for (const PolarSite& source : fg_seg) {
      for (const auto& entry : fg_sites) {
        e_erf += inter.CalcErfStaticEnergy<PolarSite, PolarSite>(
            source, *entry.first, Eigen::Vector3d::Zero());
      }
    }
    const double reference = e_real + e_recip - e_erf;

    // ---- oracle: phi and E sampled with unit probes.
    //
    // The segment id handed to AddFieldAt is a FOREGROUND one on purpose.
    // It only decides whether the zero-translation self-pair is skipped,
    // and that skip is disabled for sources that have a foreground copy
    // -- which is the behaviour a point that is not a site of its own
    // wants, and what makes this reproduce the reference.
    std::vector<std::unique_ptr<PolarSite>> unit_probes;
    unit_probes.reserve(std::size_t(n_sites) + 1);
    auto sample = [&](const Eigen::Vector3d& point, Eigen::Vector3d& field) {
      unit_probes.push_back(std::make_unique<PolarSite>(0, "H", point));
      PolarSite& probe = *unit_probes.back();
      probe.Reset();
      probe.setCharge(1.0);
      probe.setStaticDipole(Eigen::Vector3d::Zero());
      probe.setInduced_Dipole(Eigen::Vector3d::Zero());

      real_sum.AddFieldAt<Estatic::V>(fg_id, probe, EwaldChargeState::Neutral);
      double phi =
          real_sum.CalcStaticEnergyAt(probe, EwaldChargeState::Neutral);

      recip_sum.AddFieldAt<Estatic::V>(probe, EwaldChargeState::Neutral);
      const std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> one{
          {&probe, point}};
      phi += recip_sum.CalcStaticEnergyBetween(one, no_exclusions,
                                               EwaldChargeState::Neutral);

      // Subtracted for phi; ApplyErfStaticFieldCorrection already carries
      // the minus for the field, so V() ends up holding the same
      // combination.
      for (const PolarSite& source : fg_seg) {
        phi -= inter.CalcErfStaticEnergy<PolarSite, PolarSite>(
            source, probe, Eigen::Vector3d::Zero());
        inter.ApplyErfStaticFieldCorrection<PolarSite, Estatic::V>(
            source, probe, Eigen::Vector3d::Zero());
      }
      field = probe.V();
      return phi;
    };

    double oracle = 0.0;
    for (const PolarSite& site : fg_seg) {
      Eigen::Vector3d field;
      const double phi = sample(site.getPos(), field);
      oracle += site.getCharge() * phi - site.getStaticDipole().dot(field);
    }

    BOOST_REQUIRE_GT(std::abs(reference), 1e-10);
    BOOST_CHECK_SMALL(std::abs(oracle - reference) / std::abs(reference),
                      1e-10);

    Eigen::Vector3d far_field;
    const double phi_far = sample(far_point, far_field);
    if (!have_far) {
      phi_far_reference = phi_far;
      have_far = true;
      BOOST_REQUIRE_GT(std::abs(phi_far), 1e-10);
    } else {
      BOOST_CHECK_SMALL(
          std::abs(phi_far - phi_far_reference) / std::abs(phi_far_reference),
          1e-8);
    }
  }
}

// The batched potential must agree, point for point, with what a unit
// test charge reports -- the oracle the case above pinned to the energy
// convention. Same numbers by two routes: one rebuilding S(k) per point,
// one paying for it once.
//
// AGAINST THE SUM OF BOTH CHANNELS. TotalStructureFactors builds its
// moments from getStaticDipole() + getInducedDipole(), so PotentialAtMany
// is the counterpart of AddFieldAtMany and already carries the permanent
// and induced backgrounds together. Comparing it against
// CalcStaticEnergyBetween alone is wrong by the induced term, and the
// lattice below is built so that error would be tens of percent rather
// than something a loose tolerance could absorb.
BOOST_AUTO_TEST_CASE(batched_potential_matches_the_unit_probe) {
  const double L = 20.0;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;
  const Index n = 3;
  const double d = L / double(n);
  Index id = 0;
  for (Index a = 0; a < n; ++a) {
    for (Index b = 0; b < n; ++b) {
      for (Index c = 0; c < n; ++c) {
        const Eigen::Vector3d centre(double(a) * d, double(b) * d,
                                     double(c) * d);
        PolarSegment seg("seg", id);
        const double t = 0.63;
        const Eigen::Vector3d offsets[3] = {
            {0.0, 0.0, 0.0}, {t, t, t}, {-t, -t, t}};
        const double charges[3] = {-0.4, 0.2, 0.2};
        for (Index j = 0; j < 3; ++j) {
          PolarSite site(j, (j == 0) ? "C" : "H", centre + offsets[j]);
          site.setpolarization(3.0 * Eigen::Matrix3d::Identity());
          site.setCharge(charges[j]);
          // Static AND induced, neither symmetric across the cell, so
          // neither channel can pass by cancelling to zero.
          site.setStaticDipole(Eigen::Vector3d(
              0.15 * double(j + 1) + 1e-2 * double(a), -0.10 + 1e-2 * double(b),
              0.05 * double(j) - 1e-2 * double(c)));
          site.setInduced_Dipole(Eigen::Vector3d(
              1e-2 * double(j + 1) + 1e-3 * double(a), -7e-3 + 2e-3 * double(b),
              4e-3 * double(j) - 1e-3 * double(c)));
          seg.push_back(site);
        }
        registry.Register(id, EwaldChargeState::Neutral, seg);
        ++id;
      }
    }
  }

  // On a site, between segments, and off-lattice.
  const std::vector<Eigen::Vector3d> points = {{0.0, 0.0, 0.0},
                                               {0.63, 0.63, 0.63},
                                               {0.5 * d, 0.37 * d, 0.71 * d},
                                               {1.5 * d, 0.5 * d, 2.5 * d},
                                               {-0.4 * d, 1.2 * d, 0.1 * d}};

  const double alpha = 0.25;
  const double k_max = 14.0 * alpha;
  EwaldReciprocalSpaceSum recip(box, registry, alpha, k_max);

  const Eigen::VectorXd batched =
      recip.PotentialAtMany(points, EwaldChargeState::Neutral);
  BOOST_REQUIRE_EQUAL(batched.size(), Index(points.size()));

  const std::vector<const PolarSite*> no_exclusions;
  for (std::size_t p = 0; p < points.size(); ++p) {
    PolarSite probe(0, "H", points[p]);
    probe.Reset();
    probe.setCharge(1.0);
    probe.setStaticDipole(Eigen::Vector3d::Zero());
    probe.setInduced_Dipole(Eigen::Vector3d::Zero());
    const std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> one{
        {&probe, points[p]}};

    const double oracle_static = recip.CalcStaticEnergyBetween(
        one, no_exclusions, EwaldChargeState::Neutral);
    const double oracle_induced = recip.CalcInducedSourceEnergyBetween(
        one, no_exclusions, EwaldChargeState::Neutral);
    const double oracle = oracle_static + oracle_induced;

    // Not vacuous, and the induced half must be big enough that leaving
    // it out would fail this case rather than slip under the tolerance.
    BOOST_REQUIRE_GT(std::abs(oracle), 1e-12);
    BOOST_REQUIRE_GT(std::abs(oracle_induced / oracle), 1e-6);

    BOOST_CHECK_SMALL(std::abs(batched[Index(p)] - oracle) / std::abs(oracle),
                      1e-10);
  }
}

BOOST_AUTO_TEST_SUITE_END()
