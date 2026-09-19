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

#define BOOST_TEST_MODULE ewaldrealspacesum_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/ewaldrealspaceinteractor.h"
#include "votca/xtp/ewaldrealspacesum.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(ewaldrealspacesum_test)

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

// With a large box and a moderate alpha, erfc(alpha*r) underflows to
// exactly 0 in double precision for any translation other than (0,0,0),
// so the only translation that can contribute anything is the untranslated
// one. This isolates AddFieldAt's wiring (registry lookup, intermolecular
// scope, shell loop termination) from the periodic-image summation itself,
// and cross-checks it against a single direct call to the already-
// validated EwaldRealSpaceInteractor.
BOOST_AUTO_TEST_CASE(single_relevant_image_matches_direct_interactor) {
  const double alpha = 0.3;
  Eigen::Matrix3d box = 1000.0 * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;
  registry.Register(1, EwaldChargeState::Neutral,
                    MakeChargeSegment(1, Eigen::Vector3d::Zero(), 0.7));

  PolarSite target(2, "H", Eigen::Vector3d(3.0, 0.0, 0.0));

  EwaldRealSpaceSum sum(box, registry, alpha, 0.39, /*r_min=*/1.0,
                        /*field_tol=*/1e-14);
  sum.AddFieldAt<Estatic::V>(2, target, EwaldChargeState::Neutral);

  PolarSite target_direct(3, "H", Eigen::Vector3d(3.0, 0.0, 0.0));
  EwaldRealSpaceInteractor interactor(alpha, 0.39);
  PolarSite source(1, "H", Eigen::Vector3d::Zero());
  Vector9d mpoles = Vector9d::Zero();
  mpoles(0) = 0.7;
  source.setMultipole(mpoles, 0);
  interactor.ApplyStaticField<PolarSite, Estatic::V>(source, target_direct);

  BOOST_CHECK(target.V().isApprox(target_direct.V(), 1e-12));
  BOOST_CHECK(target.V().norm() > 1e-8);  // sanity: not vacuously zero
}

// A target within its own registered segment must never feel that
// segment's own field, at any translation -- this is the intermolecular-
// only scope decision documented on the class itself.
BOOST_AUTO_TEST_CASE(intramolecular_contribution_excluded) {
  const double alpha = 0.3;
  Eigen::Matrix3d box = 1000.0 * Eigen::Matrix3d::Identity();

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

  // target IS site_b, registered under the same segment id (1) it is
  // itself part of.
  PolarSite target = site_b;
  target.Reset();

  EwaldRealSpaceSum sum(box, registry, alpha, 0.39, /*r_min=*/1.0,
                        /*field_tol=*/1e-14);
  sum.AddFieldAt<Estatic::V>(1, target, EwaldChargeState::Neutral);

  BOOST_CHECK(target.V().isApprox(Eigen::Vector3d::Zero(), 1e-14));
}

// A genuinely small box, so multiple periodic images fall within
// convergence range. Cross-checks AddFieldAt's shell/convergence loop
// against a manual sum over the same explicit translations, computed with
// direct calls to EwaldRealSpaceInteractor -- this specifically exercises
// the shell-grouping and multi-shell accumulation logic that the previous
// test (single relevant image) does not.
BOOST_AUTO_TEST_CASE(periodic_images_match_manual_sum) {
  const double alpha = 0.5;
  const double L = 4.0;
  Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;
  registry.Register(1, EwaldChargeState::Neutral,
                    MakeChargeSegment(1, Eigen::Vector3d::Zero(), 0.4));

  PolarSite target(2, "H", Eigen::Vector3d(1.5, 0.3, -0.2));

  EwaldRealSpaceSum sum(box, registry, alpha, 0.39, /*r_min=*/6.0,
                        /*field_tol=*/1e-12, /*shell_width=*/0.5,
                        /*n_max=*/10);
  sum.AddFieldAt<Estatic::V>(2, target, EwaldChargeState::Neutral);

  // Manual reference: sum over every translation with |na|,|nb|,|nc|<=6
  // (comfortably enough to have converged given alpha=0.5 and the box
  // size above -- r_min=6 bohr already forces AddFieldAt itself past
  // this range).
  PolarSite target_ref(3, "H", Eigen::Vector3d(1.5, 0.3, -0.2));
  EwaldRealSpaceInteractor interactor(alpha, 0.39);
  Vector9d mpoles = Vector9d::Zero();
  mpoles(0) = 0.4;
  for (Index na = -6; na <= 6; ++na) {
    for (Index nb = -6; nb <= 6; ++nb) {
      for (Index nc = -6; nc <= 6; ++nc) {
        Eigen::Vector3d t(double(na) * L, double(nb) * L, double(nc) * L);
        PolarSite source(1, "H", t);
        source.setMultipole(mpoles, 0);
        interactor.ApplyStaticField<PolarSite, Estatic::V>(source,
                                                            target_ref);
      }
    }
  }

  BOOST_CHECK(target.V().isApprox(target_ref.V(), 1e-8));
  if (!target.V().isApprox(target_ref.V(), 1e-8)) {
    std::cout << "AddFieldAt: " << target.V().transpose() << std::endl;
    std::cout << "manual sum: " << target_ref.V().transpose() << std::endl;
  }
}


// Declaring a foreground copy must remove EXACTLY that copy's
// contribution from the periodic background -- the "carving out" an
// MM/MM or QM/MM job performs before handling those segments explicitly.
//
// Asserted as an identity against a directly computed pair sum, not
// against stored numbers: the removed amount is independently
// recomputable, so a reference cannot agree with a mistake in the
// exclusion logic. The entry count is checked too, because a correct
// TOTAL could still hide removing one image while adding another.
BOOST_AUTO_TEST_CASE(foreground_suppression_removes_exactly_one_copy) {
  const double alpha = 0.3;
  const double thole = 0.39;
  const double L = 14.0;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;
  const Index n = 2;
  const double d = L / double(n);
  Index id = 0;
  for (Index a = 0; a < n; ++a) {
    for (Index b = 0; b < n; ++b) {
      for (Index c = 0; c < n; ++c) {
        const Eigen::Vector3d centre(double(a) * d, double(b) * d,
                                     double(c) * d);
        PolarSegment seg("seg", id);
        const double t = 0.63;
        const Eigen::Vector3d offsets[5] = {{0.0, 0.0, 0.0},
                                            {t, t, t},
                                            {t, -t, -t},
                                            {-t, t, -t},
                                            {-t, -t, t}};
        for (Index j = 0; j < 5; ++j) {
          PolarSite site(j, (j == 0) ? "C" : "H", centre + offsets[j]);
          site.setpolarization(((j == 0) ? 8.0 : 3.0) *
                               Eigen::Matrix3d::Identity());
          site.setCharge((j == 0) ? -0.4 : 0.1);
          site.setStaticDipole(
              Eigen::Vector3d(1e-2 * double(j + 1), -5e-3, 2e-3 * double(j)));
          site.setInduced_Dipole(Eigen::Vector3d(1e-3, -5e-4, 7e-4) *
                                 (1.0 + 0.1 * double(j)));
          seg.push_back(site);
        }
        registry.Register(id, EwaldChargeState::Neutral, seg);
        ++id;
      }
    }
  }

  const Index target_id = 0;
  const Index fg_id = 3;
  PolarSegment& fg_seg = registry.Get(fg_id, EwaldChargeState::Neutral);
  Eigen::Vector3d fg_pos = Eigen::Vector3d::Zero();
  Index n_sites = 0;
  for (const PolarSite& site : fg_seg) {
    fg_pos += site.getPos();
    ++n_sites;
  }
  fg_pos /= double(n_sites);

  EwaldRealSpaceSum plain(box, registry, alpha, thole, 12.0, 1e-12);
  PolarSite t_plain = registry.Get(target_id, EwaldChargeState::Neutral)[0];
  t_plain.Reset();
  plain.AddFieldAt<Estatic::V>(target_id, t_plain, EwaldChargeState::Neutral);

  std::vector<std::pair<Index, Eigen::Vector3d>> foreground{{fg_id, fg_pos}};
  EwaldRealSpaceSum carved(box, registry, alpha, thole, 12.0, 1e-12, 0.945, 15,
                           6.0, foreground);
  PolarSite t_carved = registry.Get(target_id, EwaldChargeState::Neutral)[0];
  t_carved.Reset();
  carved.AddFieldAt<Estatic::V>(target_id, t_carved,
                                EwaldChargeState::Neutral);

  // The contribution of exactly that copy, computed directly.
  EwaldRealSpaceInteractor interactor(alpha, thole);
  PolarSite t_direct = registry.Get(target_id, EwaldChargeState::Neutral)[0];
  t_direct.Reset();
  for (const PolarSite& source : fg_seg) {
    interactor.ApplyStaticField<PolarSite, Estatic::V>(source, t_direct);
    interactor.ApplyInducedField<Estatic::V>(source, t_direct);
  }

  const Eigen::Vector3d removed = t_plain.V() - t_carved.V();
  BOOST_CHECK_SMALL((removed - t_direct.V()).norm() / t_direct.V().norm(),
                    1e-12);

  const auto s_plain = plain.GetNeighborStats();
  const auto s_carved = carved.GetNeighborStats();
  BOOST_CHECK_GT(s_carved.foreground, 0);
  // One copy suppressed, and nothing else gained or lost with it.
  BOOST_CHECK_EQUAL(s_plain.entries - s_carved.entries, s_carved.foreground);
}


// The accumulated Q-Q energy must use exactly the neighbour set the
// field used: same foreground suppression, same distance cull. Checked
// by suppressing one copy and requiring the energy to drop by precisely
// that copy's erfc-screened energy -- independently recomputed through
// the interactor, so a reference cannot agree with a mistake in the
// accumulation.
//
// Note the energy drops by the ERFC part only, not the bare part. That
// is correct and worth stating: suppression removes the copy from the
// real-space sum, while its reciprocal-space (erf) share is removed
// separately by the erf energy correction, exactly as for the field.
BOOST_AUTO_TEST_CASE(static_energy_uses_the_same_suppressed_neighbour_set) {
  const double alpha = 0.3;
  const double thole = 0.39;
  const double L = 14.0;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;
  const Index n = 2;
  const double d = L / double(n);
  Index id = 0;
  for (Index a = 0; a < n; ++a) {
    for (Index b = 0; b < n; ++b) {
      for (Index c = 0; c < n; ++c) {
        const Eigen::Vector3d centre(double(a) * d, double(b) * d,
                                     double(c) * d);
        PolarSegment seg("seg", id);
        const double t = 0.63;
        const Eigen::Vector3d offsets[5] = {{0.0, 0.0, 0.0},
                                            {t, t, t},
                                            {t, -t, -t},
                                            {-t, t, -t},
                                            {-t, -t, t}};
        for (Index j = 0; j < 5; ++j) {
          PolarSite site(j, (j == 0) ? "C" : "H", centre + offsets[j]);
          site.setpolarization(((j == 0) ? 8.0 : 3.0) *
                               Eigen::Matrix3d::Identity());
          site.setCharge((j == 0) ? -0.4 : 0.1);
          site.setStaticDipole(
              Eigen::Vector3d(1e-2 * double(j + 1), -5e-3, 2e-3 * double(j)));
          site.setInduced_Dipole(Eigen::Vector3d(1e-3, -5e-4, 7e-4));
          seg.push_back(site);
        }
        registry.Register(id, EwaldChargeState::Neutral, seg);
        ++id;
      }
    }
  }

  const Index target_id = 0;
  const Index fg_id = 3;
  const PolarSegment& fg_seg = registry.Get(fg_id, EwaldChargeState::Neutral);
  Eigen::Vector3d fg_pos = Eigen::Vector3d::Zero();
  Index n_sites = 0;
  for (const PolarSite& site : fg_seg) {
    fg_pos += site.getPos();
    ++n_sites;
  }
  fg_pos /= double(n_sites);

  EwaldRealSpaceSum plain(box, registry, alpha, thole, 12.0, 1e-12);
  PolarSite t_plain = registry.Get(target_id, EwaldChargeState::Neutral)[0];
  t_plain.Reset();
  plain.AddFieldAt<Estatic::V>(target_id, t_plain, EwaldChargeState::Neutral);
  const double e_plain =
      plain.CalcStaticEnergyAt(t_plain, EwaldChargeState::Neutral);

  std::vector<std::pair<Index, Eigen::Vector3d>> foreground{{fg_id, fg_pos}};
  EwaldRealSpaceSum carved(box, registry, alpha, thole, 12.0, 1e-12, 0.945, 15,
                           6.0, foreground);
  PolarSite t_carved = registry.Get(target_id, EwaldChargeState::Neutral)[0];
  t_carved.Reset();
  carved.AddFieldAt<Estatic::V>(target_id, t_carved,
                                EwaldChargeState::Neutral);
  const double e_carved =
      carved.CalcStaticEnergyAt(t_carved, EwaldChargeState::Neutral);

  EwaldRealSpaceInteractor interactor(alpha, thole);
  double e_copy = 0.0;
  for (const PolarSite& source : fg_seg) {
    e_copy +=
        interactor.CalcStaticEnergy<PolarSite, PolarSite>(source, t_carved);
  }

  BOOST_CHECK_SMALL(std::abs((e_plain - e_carved) - e_copy) / std::abs(e_copy),
                    1e-12);
}

// Querying the energy without a neighbour list must fail loudly rather
// than silently building one: that would turn a cheap query into the
// full shell search, and make the cost depend on call order.
BOOST_AUTO_TEST_CASE(static_energy_requires_an_existing_neighbour_list) {
  const double L = 14.0;
  const Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry;
  PolarSegment seg("seg", 0);
  PolarSite site(0, "C", Eigen::Vector3d::Zero());
  site.setpolarization(Eigen::Matrix3d::Identity());
  site.setCharge(-0.4);
  seg.push_back(site);
  registry.Register(0, EwaldChargeState::Neutral, seg);

  EwaldRealSpaceSum sum(box, registry, 0.3, 0.39, 12.0, 1e-12);
  PolarSite probe = registry.Get(0, EwaldChargeState::Neutral)[0];
  BOOST_CHECK_THROW(sum.CalcStaticEnergyAt(probe, EwaldChargeState::Neutral),
                    std::exception);
}

BOOST_AUTO_TEST_SUITE_END()
