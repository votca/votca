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

BOOST_AUTO_TEST_SUITE_END()
