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

// The combined real- plus reciprocal-space field (at well-converged
// cutoffs) must not depend on the Ewald splitting parameter alpha. The
// shape/k=0 term is alpha-independent and not needed here. This
// sign-sensitive check is not possible for the real-space sum alone.
BOOST_AUTO_TEST_CASE(alpha_independence_of_combined_field) {
  Eigen::Matrix3d box = 30.0 * Eigen::Matrix3d::Identity();

  auto compute_combined_field = [&](double alpha) {
    EwaldRegistry registry;
    registry.Register(1, EwaldChargeState::Neutral,
                      MakeChargeSegment(1, Eigen::Vector3d::Zero(), 0.6));

    PolarSite target(2, "H", Eigen::Vector3d(2.5, 1.0, -0.7));

    // k_max = 12*alpha keeps the truncation error fixed at
    // exp(-(12 alpha)^2/(4 alpha^2)) = exp(-36) = 2.3e-16 for every
    // alpha, so it cannot mimic an alpha dependence. It is also much
    // cheaper than a fixed k_max = 25 (98,384 vs 7,123,810 k-vectors at
    // alpha = 0.5 in this 30 bohr box).
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
// source charge along that axis, matching the EwaldRealSpaceInteractor
// sign convention.
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

// A site may be listed more than once in a batch; every occurrence must
// be added, as if AddFieldAt were called once per entry.
BOOST_AUTO_TEST_CASE(add_field_at_many_handles_repeated_targets) {
  Eigen::Matrix3d box = 10.0 * Eigen::Matrix3d::Identity();
  EwaldRegistry registry;
  registry.Register(1, EwaldChargeState::Neutral,
                    MakeChargeSegment(1, Eigen::Vector3d::Zero(), 1.0));
  EwaldReciprocalSpaceSum recip(box, registry, 0.35, /*k_max=*/10.0);

  const Index n_sites = 40;
  const Index copies = 4;
  std::vector<PolarSite> serial_sites;
  std::vector<PolarSite> batched_sites;
  for (Index i = 0; i < n_sites; ++i) {
    const Eigen::Vector3d pos(0.2 * double(i % 7) + 1.0,
                              0.3 * double(i % 5) - 0.5, 0.1 * double(i));
    serial_sites.emplace_back(10 + i, "H", pos);
    batched_sites.emplace_back(10 + i, "H", pos);
  }
  std::vector<PolarSite*> batch;
  for (Index c = 0; c < copies; ++c) {
    for (Index i = 0; i < n_sites; ++i) {
      batch.push_back(&batched_sites[std::size_t(i)]);
      recip.AddFieldAt<Estatic::V>(serial_sites[std::size_t(i)],
                                   EwaldChargeState::Neutral);
    }
  }
  recip.AddFieldAtMany<Estatic::V>(batch, EwaldChargeState::Neutral);
  for (Index i = 0; i < n_sites; ++i) {
    BOOST_CHECK(serial_sites[std::size_t(i)].V().isApprox(
        batched_sites[std::size_t(i)].V(), 1e-12));
  }
}

// Unlike EwaldRealSpaceSum, EwaldReciprocalSpaceSum excludes nothing from
// the structure factor, not even the target's own segment, as in legacy
// PolarBackground (SP_SFactorCalc/FP_KFieldCalc and
// SU_SFactorCalc/FU_KFieldCalc). The erf part of pairs excluded from the
// real-space sum must be present for erfc + erf = 1 to hold.
//
// Checks that a target feels the reciprocal-space field of another site
// of its own segment. A point charge's own contribution to its field is
// zero by construction: its structure-factor term evaluated at its own
// position is real, and the field is taken from the imaginary part.
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

  // Nonzero: site_b feels site_a (opposite charges 2.0 bohr apart).
  BOOST_CHECK(target.V().norm() > 1e-3);

  // Cross-check against a registry containing only site_a: since site_b's
  // self-contribution is zero, the fields must agree exactly.
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
// This is far more sensitive than checking either half alone. For
// example, subtracting an erf correction from the energy by analogy with
// the field (where the reciprocal part contains the total structure
// factor) would remove a foreground-foreground term that the energy's
// reciprocal cross term never contains; being erf-weighted, that error
// is alpha-dependent and fails here.
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
    // Nothing is held out of S_bg. Real space suppresses the foreground
    // segment's COINCIDENT copy, which is removed from the reciprocal
    // side by subtracting its erf energy; its other periodic images are
    // ordinary background and belong in both sums. Holding the segment
    // out of S_bg would remove those images from the reciprocal side only
    // and make the total alpha-dependent. See EwaldRegion's step (5).
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
// checks the real/reciprocal split rather than the damping
// model. A second case below measures how large the damped
// alpha-dependence actually is, so the approximation is quantified
// rather than assumed small.
BOOST_AUTO_TEST_CASE(induced_source_energy_is_independent_of_the_splitting) {
  const double L = 20.0;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  // Large enough that au3 >= 40 for every relevant pair (no damping).
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
// split is only approximately alpha-independent. At realistic damping
// the difference is negligible: damping acts only at sub-molecular
// separations, and the real-space sum excludes same-segment pairs. Measured on
// this fixture: thole_a = 0.39 gives the same energy as thole_a = 1e6 to nine
// digits, with an alpha drift of 7e-10 against 1.6e-14 for the strictly
// undamped case.
//
// The first case below checks this; the second checks that damping is
// applied at all, by making it strong enough to matter.
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

// Strong damping must change the term: at thole_a = 0.01 it changes
// sign on this fixture.
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

// A foreground that is both charged and multi-segment. Only this
// combination detects holding the whole foreground out of S_bg, which
// deletes its periodic images (a lattice of vacancies instead of one
// cavity; 1e-3 eV alpha-dependence on a production job).
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
    // Down to alpha = 0.12: for alpha >= 0.25 in this box the
    // erfc-screened part of a segment's own-image interaction is
    // negligible, and the real-space half would not be tested.
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
        // 1e-5 rather than 1e-8 as in the single-segment cases. The
        // residual shrinks ~30x per alpha step, indicating truncation
        // rather than a missing term; its origin is not understood.
        BOOST_CHECK_SMALL(std::abs(total - reference) / std::abs(reference),
                          1e-5);
      }
    }
  }
}

// ORACLE FOR EwaldRegion::PotentialAt.
//
// The potential at a point is what a UNIT TEST CHARGE placed there would
// report as its energy: with q = 1 and mu = 0, q*phi - mu.E = phi. This
// fixes the gauge of phi, which carries an arbitrary additive constant
// (the k = 0 term, i.e. the uniform neutralising background, is
// omitted) that shifts a charged region by q * phi_0. Defining phi
// through the energy expressions gives it the same constant as all
// energies computed here.
//
// So this case asserts the round trip: phi and E sampled at each
// foreground site with a unit probe, recombined as sum_i [q_i phi_i -
// mu_i . E_i], must reproduce the energy assembled directly from those
// sites, to round-off.
//
// Probes must have distinct, stable addresses: EwaldRealSpaceSum keys
// its neighbour cache on the target address, so a probe reused at
// different positions would get the first point's neighbour list. Hence
// the unique_ptrs (and PotentialAt does not use probes).
BOOST_AUTO_TEST_CASE(unit_probe_potential_reproduces_the_static_energy) {
  const double L = 20.0;
  const double thole = 0.39;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();

  // The lattice of static_energy_is_independent_of_the_splitting, but
  // with dipoles of the magnitude of a rank-1 .mps, so that the dipole
  // part of phi is exercised.
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

  // Away from the foreground segment; phi must not depend on alpha here
  // either.
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
    // A FOREGROUND segment id is passed to AddFieldAt on purpose: it only
    // controls the zero-translation self-skip, which is disabled for
    // sources with a foreground copy, as needed for a free point.
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

// The batched potential must agree point by point with a unit test
// charge (the oracle of the case above): S(k) rebuilt per point versus
// computed once.
//
// Compared against the sum of both channels: TotalStructureFactors uses
// getStaticDipole() + getInducedDipole(), so PotentialAtMany contains
// the permanent and induced backgrounds. The lattice is chosen so that
// omitting the induced term would be a tens-of-percent error.
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

    // The induced part must be large enough that omitting it would fail.
    BOOST_REQUIRE_GT(std::abs(oracle), 1e-12);
    BOOST_REQUIRE_GT(std::abs(oracle_induced / oracle), 1e-6);

    BOOST_CHECK_SMALL(std::abs(batched[Index(p)] - oracle) / std::abs(oracle),
                      1e-10);
  }
}

BOOST_AUTO_TEST_SUITE_END()
