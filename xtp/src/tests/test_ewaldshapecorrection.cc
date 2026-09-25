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
  double expected_z = -(4.0 * kPi / volume) * ExpectedTotalDipoleMoment().z();

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
// ---------------------------------------------------------------------
// Shape/surface ENERGY (the cross term), CalcStaticEnergyBetween.
// ---------------------------------------------------------------------

namespace {

constexpr double kPiLocal = 3.14159265358979323846;

PolarSite MakeSite(Index id, const Eigen::Vector3d& pos, double charge,
                   const Eigen::Vector3d& static_dipole) {
  PolarSite site(id, "H", pos);
  site.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles = Vector9d::Zero();
  mpoles(0) = charge;
  // PolarSite's Q11c, Q11s, Q10 -> x, y, z, as the registry builder
  // above already relies on.
  mpoles(1) = static_dipole.x();
  mpoles(2) = static_dipole.y();
  mpoles(3) = static_dipole.z();
  site.setMultipole(mpoles, 1);
  return site;
}

// A background with NO induced dipoles, so the energy (permanent only)
// and the field (permanent + induced) describe the same moments and can
// be checked directly against each other.
EwaldRegistry BuildPermanentOnlyBackground() {
  EwaldRegistry registry;

  PolarSegment seg1("seg", 1);
  seg1.push_back(MakeSite(1, Eigen::Vector3d(1.0, 0.0, 0.0), 0.4,
                          Eigen::Vector3d(0.0, 0.2, 0.0)));
  registry.Register(1, EwaldChargeState::Neutral, seg1);

  PolarSegment seg2("seg", 2);
  seg2.push_back(MakeSite(2, Eigen::Vector3d(0.0, 2.0, -1.0), -0.4,
                          Eigen::Vector3d(0.0, 0.0, 0.3)));
  registry.Register(2, EwaldChargeState::Neutral, seg2);

  return registry;
}

// The same three moments the implementation builds, computed here by an
// independent route so the formula tests do not simply re-run it.
struct RefMoments {
  double q0 = 0.0;
  Eigen::Vector3d q1 = Eigen::Vector3d::Zero();
  Eigen::Matrix3d q2 = Eigen::Matrix3d::Zero();
};

RefMoments MomentsOf(const PolarSegment& seg, const Eigen::Vector3d& shift) {
  RefMoments m;
  for (const PolarSite& site : seg) {
    const Eigen::Vector3d r = site.getPos() + shift;
    const double q = site.getCharge();
    const Eigen::Vector3d mu = site.getStaticDipole();
    m.q0 += q;
    m.q1 += q * r + mu;
    m.q2 += 0.5 * q * r * r.transpose() + mu * r.transpose();
  }
  return m;
}

std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> AsForeground(
    const PolarSegment& seg,
    const Eigen::Vector3d& shift = Eigen::Vector3d::Zero()) {
  std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> sites;
  for (const PolarSite& site : seg) {
    sites.push_back({&site, site.getPos() + shift});
  }
  return sites;
}

// A charged foreground carrying a dipole as well, so every term of the
// bracket contributes.
PolarSegment ChargedForeground() {
  PolarSegment fg("fg", 7);
  fg.push_back(MakeSite(10, Eigen::Vector3d(3.0, 1.0, -2.0), 1.0,
                        Eigen::Vector3d(0.15, -0.1, 0.05)));
  fg.push_back(MakeSite(11, Eigen::Vector3d(3.5, 1.2, -2.4), -0.3,
                        Eigen::Vector3d(0.0, 0.25, 0.0)));
  return fg;
}

}  // namespace

// THE test for this term. The shape energy must not depend on where the
// origin is, and for a CHARGED foreground that is a real constraint: the
// dipole moments of both sides shift under translation, and only the
// full three-term bracket cancels those shifts.
//
// The Q1.Q1-only form this method first had passes every other test in
// this file and fails this one.
BOOST_AUTO_TEST_CASE(shape_energy_is_origin_invariant) {
  const double volume = 125.0;
  const PolarSegment fg = ChargedForeground();

  // An arbitrary, deliberately un-symmetric shift, comparable in size to
  // the coordinates themselves so any residual dependence is visible.
  const Eigen::Vector3d a(7.3, -4.1, 2.9);

  for (EwaldShape shape_kind : {EwaldShape::Cube, EwaldShape::Slab}) {
    EwaldRegistry unshifted = BuildPermanentOnlyBackground();
    EwaldShapeCorrection shape_unshifted(volume, unshifted, shape_kind);
    const double e_unshifted = shape_unshifted.CalcStaticEnergyBetween(
        AsForeground(fg), {}, EwaldChargeState::Neutral);

    // Shift EVERYTHING: the background sites are rebuilt at shifted
    // positions, and the foreground is handed the shifted positions.
    EwaldRegistry shifted;
    {
      PolarSegment seg1("seg", 1);
      seg1.push_back(MakeSite(1, Eigen::Vector3d(1.0, 0.0, 0.0) + a, 0.4,
                              Eigen::Vector3d(0.0, 0.2, 0.0)));
      shifted.Register(1, EwaldChargeState::Neutral, seg1);

      PolarSegment seg2("seg", 2);
      seg2.push_back(MakeSite(2, Eigen::Vector3d(0.0, 2.0, -1.0) + a, -0.4,
                              Eigen::Vector3d(0.0, 0.0, 0.3)));
      shifted.Register(2, EwaldChargeState::Neutral, seg2);
    }
    EwaldShapeCorrection shape_shifted(volume, shifted, shape_kind);
    const double e_shifted = shape_shifted.CalcStaticEnergyBetween(
        AsForeground(fg, a), {}, EwaldChargeState::Neutral);

    BOOST_CHECK_CLOSE(e_unshifted, e_shifted, 1e-9);
    BOOST_CHECK(std::abs(e_unshifted) > 1e-8);  // not vacuously zero
  }
}

// The full bracket, against a hand-built reference:
//
//   E = -(4*pi/3V) [ Q0_1 TrQ2_2 + Q0_2 TrQ2_1 - Q1_1 . Q1_2 ]
BOOST_AUTO_TEST_CASE(cube_shape_energy_matches_bracket_form) {
  EwaldRegistry registry = BuildPermanentOnlyBackground();
  const double volume = 125.0;
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

  const PolarSegment fg = ChargedForeground();
  const RefMoments m_fg = MomentsOf(fg, Eigen::Vector3d::Zero());

  RefMoments m_bg;
  for (Index id : {Index(1), Index(2)}) {
    const RefMoments part = MomentsOf(
        registry.Get(id, EwaldChargeState::Neutral), Eigen::Vector3d::Zero());
    m_bg.q0 += part.q0;
    m_bg.q1 += part.q1;
    m_bg.q2 += part.q2;
  }

  const double bracket = m_fg.q0 * m_bg.q2.trace() + m_bg.q0 * m_fg.q2.trace() -
                         m_fg.q1.dot(m_bg.q1);
  const double expected = -(4.0 * kPiLocal / (3.0 * volume)) * bracket;

  BOOST_CHECK_CLOSE(shape.CalcStaticEnergyBetween(AsForeground(fg), {},
                                                  EwaldChargeState::Neutral),
                    expected, 1e-9);

  // The second-moment terms must actually carry weight here, or the test
  // above would not distinguish the two forms.
  const double dipole_only =
      -(4.0 * kPiLocal / (3.0 * volume)) * (-m_fg.q1.dot(m_bg.q1));
  BOOST_CHECK(std::abs(expected - dipole_only) > 1e-6);
}

// Slab boundary: the zz components rather than the traces, with 4*pi/V
// rather than 4*pi/3V.
BOOST_AUTO_TEST_CASE(slab_shape_energy_matches_bracket_form) {
  EwaldRegistry registry = BuildPermanentOnlyBackground();
  const double volume = 125.0;
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Slab);

  const PolarSegment fg = ChargedForeground();
  const RefMoments m_fg = MomentsOf(fg, Eigen::Vector3d::Zero());

  RefMoments m_bg;
  for (Index id : {Index(1), Index(2)}) {
    const RefMoments part = MomentsOf(
        registry.Get(id, EwaldChargeState::Neutral), Eigen::Vector3d::Zero());
    m_bg.q0 += part.q0;
    m_bg.q1 += part.q1;
    m_bg.q2 += part.q2;
  }

  const double bracket = m_fg.q0 * m_bg.q2(2, 2) + m_bg.q0 * m_fg.q2(2, 2) -
                         m_fg.q1.z() * m_bg.q1.z();
  const double expected = -(4.0 * kPiLocal / volume) * bracket;

  BOOST_CHECK_CLOSE(shape.CalcStaticEnergyBetween(AsForeground(fg), {},
                                                  EwaldChargeState::Neutral),
                    expected, 1e-9);
  BOOST_CHECK(std::abs(expected) > 1e-8);
}

// For a NEUTRAL foreground the second-moment terms drop out (Q0_1 = 0,
// and the background is neutral too), and the shape energy reduces to
// the shape FIELD contracted with the foreground's dipole moment:
//
//   U = - E . Q1_fg
//
// That identity is what the term was first built from. It is true, but
// only here -- which is exactly why building on it alone was wrong.
BOOST_AUTO_TEST_CASE(neutral_foreground_reduces_to_field_contraction) {
  EwaldRegistry registry = BuildPermanentOnlyBackground();
  const double volume = 125.0;
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

  // Neutral overall, but carrying a dipole so the check is not trivial.
  PolarSegment fg("fg", 7);
  fg.push_back(MakeSite(10, Eigen::Vector3d(3.0, 1.0, -2.0), 0.6,
                        Eigen::Vector3d(0.15, -0.1, 0.05)));
  fg.push_back(MakeSite(11, Eigen::Vector3d(3.5, 1.2, -2.4), -0.6,
                        Eigen::Vector3d(0.0, 0.25, 0.0)));

  PolarSite probe(99, "X", Eigen::Vector3d::Zero());
  shape.AddFieldAt<Estatic::V>(probe, EwaldChargeState::Neutral);
  const double expected =
      -probe.V().dot(MomentsOf(fg, Eigen::Vector3d::Zero()).q1);

  BOOST_CHECK_CLOSE(shape.CalcStaticEnergyBetween(AsForeground(fg), {},
                                                  EwaldChargeState::Neutral),
                    expected, 1e-9);
  BOOST_CHECK(std::abs(expected) > 1e-8);
}

// The background's INDUCED dipoles belong to the field, not to this
// energy: the polar region accounts for them itself, from that same
// field. Switching one on must move the field and leave the energy
// exactly where it was.
BOOST_AUTO_TEST_CASE(shape_energy_ignores_background_induced_dipoles) {
  const double volume = 125.0;
  const PolarSegment fg = ChargedForeground();

  EwaldRegistry plain = BuildPermanentOnlyBackground();
  EwaldShapeCorrection shape_plain(volume, plain, EwaldShape::Cube);
  const double e_plain = shape_plain.CalcStaticEnergyBetween(
      AsForeground(fg), {}, EwaldChargeState::Neutral);

  // Same background, one induced dipole switched on.
  EwaldRegistry induced;
  {
    PolarSegment seg1("seg", 1);
    PolarSite s1 = MakeSite(1, Eigen::Vector3d(1.0, 0.0, 0.0), 0.4,
                            Eigen::Vector3d(0.0, 0.2, 0.0));
    s1.setInduced_Dipole(Eigen::Vector3d(0.7, -0.5, 0.9));
    seg1.push_back(s1);
    induced.Register(1, EwaldChargeState::Neutral, seg1);

    PolarSegment seg2("seg", 2);
    seg2.push_back(MakeSite(2, Eigen::Vector3d(0.0, 2.0, -1.0), -0.4,
                            Eigen::Vector3d(0.0, 0.0, 0.3)));
    induced.Register(2, EwaldChargeState::Neutral, seg2);
  }
  EwaldShapeCorrection shape_induced(volume, induced, EwaldShape::Cube);
  const double e_induced = shape_induced.CalcStaticEnergyBetween(
      AsForeground(fg), {}, EwaldChargeState::Neutral);

  BOOST_CHECK_CLOSE(e_plain, e_induced, 1e-12);

  // Guard against that passing because the induced dipole was too small
  // to matter: it does move the FIELD.
  PolarSite probe_a(98, "X", Eigen::Vector3d::Zero());
  shape_plain.AddFieldAt<Estatic::V>(probe_a, EwaldChargeState::Neutral);
  PolarSite probe_b(99, "X", Eigen::Vector3d::Zero());
  shape_induced.AddFieldAt<Estatic::V>(probe_b, EwaldChargeState::Neutral);
  BOOST_CHECK((probe_a.V() - probe_b.V()).norm() > 1e-6);
}

// Listed background sites are held out, matching the suppression the
// real- and reciprocal-space sums apply to the foreground's own copies.
BOOST_AUTO_TEST_CASE(shape_energy_excludes_listed_background_sites) {
  EwaldRegistry registry = BuildPermanentOnlyBackground();
  const double volume = 125.0;
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

  const PolarSegment fg = ChargedForeground();

  // Hold out segment 2 entirely; the background is then segment 1 alone.
  std::vector<const PolarSite*> excl;
  for (const PolarSite& site : registry.Get(2, EwaldChargeState::Neutral)) {
    excl.push_back(&site);
  }

  const double e_excluded = shape.CalcStaticEnergyBetween(
      AsForeground(fg), excl, EwaldChargeState::Neutral);

  const RefMoments m_fg = MomentsOf(fg, Eigen::Vector3d::Zero());
  const RefMoments m_bg = MomentsOf(registry.Get(1, EwaldChargeState::Neutral),
                                    Eigen::Vector3d::Zero());
  const double bracket = m_fg.q0 * m_bg.q2.trace() + m_bg.q0 * m_fg.q2.trace() -
                         m_fg.q1.dot(m_bg.q1);
  BOOST_CHECK_CLOSE(e_excluded, -(4.0 * kPiLocal / (3.0 * volume)) * bracket,
                    1e-9);

  // And excluding something genuinely changes the answer.
  const double e_full = shape.CalcStaticEnergyBetween(
      AsForeground(fg), {}, EwaldChargeState::Neutral);
  BOOST_CHECK(std::abs(e_full - e_excluded) > 1e-8);
}

// Regression guard for the bug this term was added alongside: the energy
// must follow the FOREGROUND's own moments. Were they taken from the
// background copies instead, two charge states of the same segment at
// the same position would report the same energy, and the difference
// between them -- the whole point of a site-energy calculation -- would
// silently lose this contribution.
BOOST_AUTO_TEST_CASE(shape_energy_tracks_foreground_charge_state) {
  EwaldRegistry registry = BuildPermanentOnlyBackground();
  const double volume = 125.0;
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

  const Eigen::Vector3d pos(3.0, 1.0, -2.0);
  const Eigen::Vector3d dipole(0.15, -0.1, 0.05);

  PolarSegment neutral_fg("fg", 7);
  neutral_fg.push_back(MakeSite(10, pos, 0.0, dipole));

  PolarSegment charged_fg("fg", 7);
  charged_fg.push_back(MakeSite(10, pos, 1.0, dipole));

  const double e_neutral = shape.CalcStaticEnergyBetween(
      AsForeground(neutral_fg), {}, EwaldChargeState::Neutral);
  const double e_charged = shape.CalcStaticEnergyBetween(
      AsForeground(charged_fg), {}, EwaldChargeState::Neutral);

  BOOST_CHECK(std::abs(e_charged - e_neutral) > 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
