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

#define BOOST_TEST_MODULE ewaldperiodicdipolesymmetry_test

// Standard includes
#include <random>

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/ewaldreciprocalspacesum.h"
#include "votca/xtp/ewaldrealspacesum.h"

using namespace votca::xtp;
using namespace votca;

// This is not a test of any one class in isolation; it is a check on
// whether the *combination* of EwaldRealSpaceSum and
// EwaldReciprocalSpaceSum, used as a periodic dipole-dipole interaction
// operator (given a trial set of induced dipoles, return the field that
// set produces at every site), is symmetric -- a hard requirement for
// using Eigen::ConjugateGradient on it later, and one that had not
// actually been checked anywhere before this test, only assumed from the
// underlying physics (interaction tensors are pairwise symmetric, and
// that ought to survive summing over periodic images, but "ought to" is
// not evidence).
//
// This also doubles as the first working exercise of the mechanism the
// eventual PCG-based induction solve will need: temporarily writing a
// trial induced-dipole vector into the registry, evaluating the periodic
// field it produces, and reading that back out as a plain vector.

BOOST_AUTO_TEST_SUITE(ewaldperiodicdipolesymmetry_test)

namespace {

PolarSegment MakeSegment(Index id, const Eigen::Vector3d& pos) {
  PolarSegment seg("seg", id);
  PolarSite site(id, "H", pos);
  site.setpolarization(Eigen::Matrix3d::Identity());
  seg.push_back(site);
  return seg;
}

// Sets every registered site's induced dipole from a flat 3N vector
// (segment order given by ids, 3 components per id, in that order),
// clears every site's field accumulators, applies the real+reciprocal
// periodic field from every *other* segment, and reads the result back
// out into a flat 3N vector -- i.e. this computes A*v for the periodic
// dipole-dipole operator A implied by real_sum/recip_sum.
Eigen::VectorXd Apply(EwaldRegistry& registry,
                     const std::vector<Index>& ids,
                     const EwaldRealSpaceSum& real_sum,
                     const EwaldReciprocalSpaceSum& recip_sum,
                     const Eigen::VectorXd& v) {
  // Write v into the registry's induced dipoles, and collect (id, site*)
  // pairs for the batched reciprocal-space call.
  std::vector<std::pair<Index, PolarSite*>> targets;
  for (std::size_t n = 0; n < ids.size(); ++n) {
    PolarSegment& seg = registry.Get(ids[n], EwaldChargeState::Neutral);
    PolarSite& site = seg[0];
    site.setInduced_Dipole(v.segment<3>(3 * Index(n)));
    site.Reset();
  }
  for (std::size_t n = 0; n < ids.size(); ++n) {
    PolarSegment& seg = registry.Get(ids[n], EwaldChargeState::Neutral);
    targets.push_back({ids[n], &seg[0]});
  }

  for (const auto& entry : targets) {
    real_sum.AddFieldAt<Estatic::V>(entry.first, *entry.second,
                                    EwaldChargeState::Neutral);
  }
  // EwaldReciprocalSpaceSum no longer takes a segment id (it never
  // excludes anything -- see its own class documentation), so only the
  // bare site pointers are needed here.
  std::vector<PolarSite*> recip_targets;
  recip_targets.reserve(targets.size());
  for (const auto& entry : targets) {
    recip_targets.push_back(entry.second);
  }
  recip_sum.AddFieldAtMany<Estatic::V>(recip_targets,
                                       EwaldChargeState::Neutral);

  Eigen::VectorXd result(v.size());
  for (std::size_t n = 0; n < ids.size(); ++n) {
    PolarSegment& seg = registry.Get(ids[n], EwaldChargeState::Neutral);
    result.segment<3>(3 * Index(n)) = seg[0].V();
  }
  return result;
}

}  // namespace

BOOST_AUTO_TEST_CASE(periodic_dipole_operator_is_symmetric) {
  Eigen::Matrix3d box = 12.0 * Eigen::Matrix3d::Identity();
  const std::vector<Index> ids = {1, 2, 3};

  EwaldRegistry registry;
  registry.Register(1, EwaldChargeState::Neutral,
                    MakeSegment(1, Eigen::Vector3d(0.0, 0.0, 0.0)));
  registry.Register(2, EwaldChargeState::Neutral,
                    MakeSegment(2, Eigen::Vector3d(3.0, 1.0, -0.5)));
  registry.Register(3, EwaldChargeState::Neutral,
                    MakeSegment(3, Eigen::Vector3d(-2.0, 2.5, 1.5)));

  EwaldRealSpaceSum real_sum(box, registry, /*alpha=*/0.35, /*thole_a=*/0.39,
                            /*r_min=*/6.0, /*field_tol=*/1e-12);
  EwaldReciprocalSpaceSum recip_sum(box, registry, /*alpha=*/0.35,
                                    /*k_max=*/20.0);

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  Eigen::VectorXd u(9), v(9);
  for (Index i = 0; i < 9; ++i) {
    u(i) = dist(rng);
    v(i) = dist(rng);
  }

  Eigen::VectorXd Av = Apply(registry, ids, real_sum, recip_sum, v);
  Eigen::VectorXd Au = Apply(registry, ids, real_sum, recip_sum, u);

  double uAv = u.dot(Av);
  double vAu = v.dot(Au);

  BOOST_CHECK_CLOSE(uAv, vAu, 1e-4);
  // Sanity: not vacuously symmetric because both are ~zero.
  BOOST_CHECK(std::abs(uAv) > 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
