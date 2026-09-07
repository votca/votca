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

#define BOOST_TEST_MODULE ewaldshapecorrection_test

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/ewaldshapecorrection.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(ewaldshapecorrection_test)

// A small system with a charge, a static dipole, and an induced dipole,
// spread across two registered segments, chosen so the total dipole
// moment M = sum(q*r + mu_static + mu_induced) is nonzero and has no
// special symmetry -- exercising every term of TotalDipoleMoment's own
// sum, and giving a genuine (not vacuously zero) check for the field
// formulas below.
EwaldRegistry BuildTestRegistry() {
  EwaldRegistry registry;

  PolarSegment seg1("seg", 1);
  PolarSite site1(1, "H", Eigen::Vector3d(1.0, 0.0, 0.0));
  site1.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles1 = Vector9d::Zero();
  mpoles1(0) = 0.5;  // charge
  site1.setMultipole(mpoles1, 1);
  site1.setInduced_Dipole(Eigen::Vector3d(0.1, 0.0, 0.0));
  seg1.push_back(site1);
  registry.Register(1, EwaldChargeState::Neutral, seg1);

  PolarSegment seg2("seg", 2);
  PolarSite site2(2, "H", Eigen::Vector3d(0.0, 2.0, 0.0));
  site2.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles2 = Vector9d::Zero();
  mpoles2(0) = -0.2;  // charge
  mpoles2(3) = 0.3;   // static dipole z-component (Q10, per PolarSite's
                      // own Q11c,Q11s,Q10 -> x,y,z convention)
  site2.setMultipole(mpoles2, 1);
  seg2.push_back(site2);
  registry.Register(2, EwaldChargeState::Neutral, seg2);

  return registry;
}

// Independent hand computation of M, to cross-check against
// EwaldShapeCorrection's own internal TotalDipoleMoment (not directly
// exposed, so this is verified indirectly through the cube-shape field
// formula in the next test instead of calling it directly).
Eigen::Vector3d ExpectedTotalDipoleMoment() {
  Eigen::Vector3d M = Eigen::Vector3d::Zero();
  // site1: q=0.5 at (1,0,0), static dipole zero, induced (0.1,0,0)
  M += 0.5 * Eigen::Vector3d(1.0, 0.0, 0.0);
  M += Eigen::Vector3d(0.1, 0.0, 0.0);
  // site2: q=-0.2 at (0,2,0), static dipole (0,0,0.3), induced zero
  M += -0.2 * Eigen::Vector3d(0.0, 2.0, 0.0);
  M += Eigen::Vector3d(0.0, 0.0, 0.3);
  return M;
}

BOOST_AUTO_TEST_CASE(cube_shape_matches_formula) {
  EwaldRegistry registry = BuildTestRegistry();
  const double volume = 125.0;

  PolarSite target(3, "H", Eigen::Vector3d(4.0, -1.0, 2.0));
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);
  shape.AddFieldAt<Estatic::V>(target, EwaldChargeState::Neutral);

  const double kPi = 3.14159265358979323846;
  Eigen::Vector3d expected =
      -(4.0 * kPi / (3.0 * volume)) * ExpectedTotalDipoleMoment();

  BOOST_CHECK(target.V().isApprox(expected, 1e-12));
  BOOST_CHECK(target.V().norm() > 1e-6);  // sanity: not vacuously zero
}

BOOST_AUTO_TEST_CASE(slab_shape_only_z_component_nonzero) {
  EwaldRegistry registry = BuildTestRegistry();
  const double volume = 125.0;

  PolarSite target(3, "H", Eigen::Vector3d(4.0, -1.0, 2.0));
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Slab);
  shape.AddFieldAt<Estatic::V>(target, EwaldChargeState::Neutral);

  const double kPi = 3.14159265358979323846;
  double expected_z =
      -(4.0 * kPi / volume) * ExpectedTotalDipoleMoment().z();

  BOOST_CHECK_EQUAL(target.V().x(), 0.0);
  BOOST_CHECK_EQUAL(target.V().y(), 0.0);
  BOOST_CHECK_CLOSE(target.V().z(), expected_z, 1e-9);
  BOOST_CHECK(std::abs(target.V().z()) > 1e-6);  // sanity: not vacuously
                                                  // zero
}

// The shape correction is a uniform (position-independent) field -- the
// depolarizing field inside a uniformly polarized sample does not depend
// on where within the sample you evaluate it.
BOOST_AUTO_TEST_CASE(field_is_position_independent) {
  EwaldRegistry registry = BuildTestRegistry();
  const double volume = 125.0;
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

  PolarSite target_a(3, "H", Eigen::Vector3d(0.0, 0.0, 0.0));
  shape.AddFieldAt<Estatic::V>(target_a, EwaldChargeState::Neutral);

  PolarSite target_b(4, "H", Eigen::Vector3d(10.0, -5.0, 3.0));
  shape.AddFieldAt<Estatic::V>(target_b, EwaldChargeState::Neutral);

  BOOST_CHECK(target_a.V().isApprox(target_b.V(), 1e-12));
}

BOOST_AUTO_TEST_SUITE_END()
