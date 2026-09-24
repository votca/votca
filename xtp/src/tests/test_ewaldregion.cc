/*
 *            Copyright 2009-2026 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE ewaldregion_test

// Standard includes
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <utility>
#include <vector>

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/property.h"
#include "votca/xtp/ewaldparameters.h"
#include "votca/xtp/ewaldregion.h"
#include "votca/xtp/ewaldregistry.h"
#include "votca/xtp/logger.h"

using namespace votca;
using namespace votca::xtp;

namespace {

// Written exactly the way EwaldBackground writes its own checkpoint:
// the registry at the top level, the parameters in an "ewald_parameters"
// child. If those two ever drift apart, this test fails -- which is the
// point, since EwaldRegion is the consumer of that layout and nothing
// else would notice a change until a real job aborted.
void WriteBackground(const std::string& file, const EwaldParameters& params,
                     Index n_segments) {
  EwaldRegistry registry;
  for (Index id = 0; id < n_segments; ++id) {
    PolarSegment seg("seg", id);
    for (Index j = 0; j < 3; ++j) {
      PolarSite site(j, (j == 0) ? "C" : "H",
                     Eigen::Vector3d(double(id) * 3.0, double(j) * 0.7,
                                     0.2 * double(j)));
      site.setpolarization(((j == 0) ? 8.0 : 3.0) *
                           Eigen::Matrix3d::Identity());
      site.setCharge((j == 0) ? -0.4 : 0.2);
      site.setInduced_Dipole(
          Eigen::Vector3d(1e-3, -2e-3, 5e-4) * (1.0 + 0.1 * double(id)));
      seg.push_back(site);
    }
    registry.Register(id, EwaldChargeState::Neutral, seg);
  }

  CheckpointFile cpf(file, CheckpointAccessLevel::CREATE);
  CheckpointWriter w = cpf.getWriter();
  registry.WriteToCpt(w);
  CheckpointWriter wp = w.openChild("ewald_parameters");
  params.WriteToCpt(wp);
}

EwaldParameters ReferenceParameters() {
  EwaldParameters params;
  params.alpha = 0.105835;
  params.k_max = 1.27003;
  params.r_min = 75.589;
  params.field_tol = 1e-12;
  params.thole_a = 0.39;
  params.screening_factor = 6.0;
  params.shape = EwaldShape::Slab;
  params.box = 66.1 * Eigen::Matrix3d::Identity();
  return params;
}

tools::Property RegionDefinition(const std::string& file) {
  tools::Property prop;
  tools::Property& region = prop.add("ewaldregion", "");
  region.add("checkpoint", file);
  return prop;
}

// A background of ONE point charge, so the field it produces at the
// origin has an unambiguous direction. Everything else in this file
// builds neutral multi-site segments, which is right for the load tests
// but useless for pinning a sign.
void WriteSingleCharge(const std::string& file, double charge,
                       const Eigen::Vector3d& position,
                       Index probe_id,
                       const Eigen::Vector3d& probe_position) {
  EwaldRegistry registry;
  // +q at `position` and -q mirrored through the probe. Two reasons, both
  // necessary:
  //
  //  * the cell must be NEUTRAL. The reciprocal sum omits k = 0, which
  //    is only legitimate for a neutral cell; a lone charge would make
  //    the periodic problem ill-posed and the field it produces
  //    meaningless to assert anything about.
  //  * mirroring makes the two contributions ADD at the probe rather
  //    than cancel: the field from +q points away from it, and the field
  //    from -q points towards it, which is the same direction. So the
  //    direction being tested is unambiguous and the signal is doubled.
  PolarSegment seg("chg", 0);
  PolarSite site(0, "H", position);
  site.setpolarization(1e-6 * Eigen::Matrix3d::Identity());
  site.setCharge(charge);
  seg.push_back(site);
  PolarSite counter(1, "H", 2.0 * probe_position - position);
  counter.setpolarization(1e-6 * Eigen::Matrix3d::Identity());
  counter.setCharge(-charge);
  seg.push_back(counter);
  registry.Register(0, EwaldChargeState::Neutral, seg);

  // The probe's own segment must also exist in the background: a
  // foreground is CARVED OUT of the background, so every foreground
  // segment has a counterpart there, and EwaldRegion rightly refuses a
  // foreground that does not. Given zero multipoles so it contributes
  // nothing to the field being measured -- its copy is suppressed from
  // the real-space sum and its (zero) erf share is removed, both
  // identically zero, leaving only the charge above.
  PolarSegment probe_seg("probe", probe_id);
  PolarSite probe_site(0, "C", probe_position);
  probe_site.setpolarization(8.0 * Eigen::Matrix3d::Identity());
  probe_site.setCharge(0.0);
  probe_seg.push_back(probe_site);
  registry.Register(probe_id, EwaldChargeState::Neutral, probe_seg);

  EwaldParameters params = ReferenceParameters();
  params.shape = EwaldShape::Cube;

  CheckpointFile cpf(file, CheckpointAccessLevel::CREATE);
  CheckpointWriter w = cpf.getWriter();
  registry.WriteToCpt(w);
  CheckpointWriter wp = w.openChild("ewald_parameters");
  params.WriteToCpt(wp);
}

// A foreground of one polarisable site, placed away from the background
// charge so the separation is well defined.
std::vector<PolarSegment> SingleProbe(Index id,
                                      const Eigen::Vector3d& position) {
  PolarSegment seg("probe", id);
  PolarSite site(0, "C", position);
  site.setpolarization(8.0 * Eigen::Matrix3d::Identity());
  site.setCharge(0.0);
  seg.push_back(site);
  return std::vector<PolarSegment>{seg};
}

}  // namespace

BOOST_AUTO_TEST_SUITE(ewaldregion_test)

BOOST_AUTO_TEST_CASE(loads_background_and_parameters_from_checkpoint) {
  const std::string file = "ewaldregion_test_background.hdf5";
  const EwaldParameters params = ReferenceParameters();
  WriteBackground(file, params, 4);

  Logger log;
  log.setReportLevel(Log::error);
  EwaldRegion region(7, log);
  tools::Property prop = RegionDefinition(file);
  region.Initialize(prop.get("ewaldregion"));

  BOOST_CHECK_EQUAL(region.identify(), "ewaldregion");
  BOOST_CHECK_EQUAL(region.size(), 4);

  // Each segment is -0.4 + 0.2 + 0.2 = 0, so the cell is neutral.
  BOOST_CHECK_SMALL(region.charge(), 1e-12);

  // The parameters must come back bit-for-bit: they are the guarantee
  // that a foreground is embedded with the same alpha the background was
  // converged with, so "close enough" is not good enough here.
  const EwaldParameters& read = region.Parameters();
  BOOST_CHECK_EQUAL(read.alpha, params.alpha);
  BOOST_CHECK_EQUAL(read.k_max, params.k_max);
  BOOST_CHECK_EQUAL(read.r_min, params.r_min);
  BOOST_CHECK_EQUAL(read.field_tol, params.field_tol);
  BOOST_CHECK_EQUAL(read.thole_a, params.thole_a);
  BOOST_CHECK_EQUAL(read.screening_factor, params.screening_factor);
  BOOST_CHECK(read.shape == params.shape);
  BOOST_CHECK_EQUAL((read.box - params.box).norm(), 0.0);

  // The converged induced dipoles are the whole reason the checkpoint
  // exists; losing them would leave a background that looks structurally
  // fine and is physically empty.
  const PolarSegment& seg =
      region.Registry().Get(2, EwaldChargeState::Neutral);
  const Eigen::Vector3d expected =
      Eigen::Vector3d(1e-3, -2e-3, 5e-4) * (1.0 + 0.1 * 2.0);
  BOOST_CHECK_SMALL((seg[0].getInducedDipole() - expected).norm(), 1e-14);

  std::remove(file.c_str());
}

// The frozen-background contract. These are not incidental: the
// inter-region SCF loop terminates only if every region reports
// Converged(), and nothing in a job re-polarizes this one, so its
// Interactwith* results must be exactly zero rather than small.
BOOST_AUTO_TEST_CASE(region_is_frozen_and_always_converged) {
  const std::string file = "ewaldregion_test_frozen.hdf5";
  WriteBackground(file, ReferenceParameters(), 2);

  Logger log;
  log.setReportLevel(Log::error);
  EwaldRegion region(3, log);
  tools::Property prop = RegionDefinition(file);
  region.Initialize(prop.get("ewaldregion"));

  BOOST_CHECK(region.Converged());
  // Reset must not disturb the loaded state -- there is no per-iteration
  // state to clear, and clearing the dipoles would silently empty the
  // background.
  const Eigen::Vector3d before =
      region.Registry().Get(1, EwaldChargeState::Neutral)[0]
          .getInducedDipole();
  region.Reset();
  BOOST_CHECK(region.Converged());
  const Eigen::Vector3d after =
      region.Registry().Get(1, EwaldChargeState::Neutral)[0]
          .getInducedDipole();
  BOOST_CHECK_SMALL((before - after).norm(), 1e-16);
  BOOST_CHECK_GT(before.norm(), 0.0);

  std::remove(file.c_str());
}

// A checkpoint without the parameter block is one written before
// EwaldParameters existed. Reading it must fail loudly: the alternative
// is defaulting the parameters and embedding a foreground with an alpha
// the background was never converged at, which nothing downstream could
// detect.
BOOST_AUTO_TEST_CASE(rejects_checkpoint_without_parameters) {
  const std::string file = "ewaldregion_test_noparams.hdf5";
  {
    EwaldRegistry registry;
    PolarSegment seg("seg", 0);
    PolarSite site(0, "C", Eigen::Vector3d::Zero());
    site.setpolarization(Eigen::Matrix3d::Identity());
    seg.push_back(site);
    registry.Register(0, EwaldChargeState::Neutral, seg);
    CheckpointFile cpf(file, CheckpointAccessLevel::CREATE);
    CheckpointWriter w = cpf.getWriter();
    registry.WriteToCpt(w);  // no ewald_parameters child
  }

  Logger log;
  log.setReportLevel(Log::error);
  EwaldRegion region(0, log);
  tools::Property prop = RegionDefinition(file);
  BOOST_CHECK_THROW(region.Initialize(prop.get("ewaldregion")),
                    std::exception);

  std::remove(file.c_str());
}

// Likewise a parameter block that is present but implausible. alpha = 0
// is not a usable Ewald splitting at all, and silently proceeding would
// produce a background whose real- and reciprocal-space halves do not
// describe the same system.
BOOST_AUTO_TEST_CASE(rejects_implausible_parameters) {
  const std::string file = "ewaldregion_test_badparams.hdf5";
  EwaldParameters bad = ReferenceParameters();
  bad.alpha = 0.0;
  WriteBackground(file, bad, 1);

  Logger log;
  log.setReportLevel(Log::error);
  EwaldRegion region(0, log);
  tools::Property prop = RegionDefinition(file);
  BOOST_CHECK_THROW(region.Initialize(prop.get("ewaldregion")),
                    std::exception);

  std::remove(file.c_str());
}


// THE SIGN AT THE REGION BOUNDARY.
//
// The Ewald code and the region framework store V with opposite signs:
// EwaldBackground builds its rhs as b = +V, PolarRegion builds its rhs
// as b = -(V + V_noE). ApplyFieldTo therefore has to hand over the field
// in the POLAR region's convention, and this asserts that it does.
//
// The quantity checked is -V, because that is literally what
// PolarRegion::CalcInducedDipolesViaPCG uses as its rhs. For a positive
// background charge, the physical field at the probe points AWAY from
// that charge, so -V must too.
//
// This is the check that was missing when the coupling first ran: the
// field was handed over un-negated, the polar region negated it again,
// and a neutral MM/MM job came out 100-700% off with a residual that was
// flat in distance rather than growing outward.
//
// CAVEAT worth knowing: the convention being asserted lives in
// PolarRegion, not here. If PolarRegion's rhs sign is ever changed, this
// test keeps passing while the coupling silently breaks again. Tying the
// two together would need a real PolarRegion in this test, which drags
// in the whole region stack.
BOOST_AUTO_TEST_CASE(field_is_handed_over_in_the_polar_region_convention) {
  const std::string file = "ewaldregion_test_sign.hdf5";
  const Eigen::Vector3d charge_pos(10.0, 0.0, 0.0);
  const Eigen::Vector3d probe_pos = Eigen::Vector3d::Zero();
  WriteSingleCharge(file, +1.0, charge_pos, 5, probe_pos);

  Logger log;
  log.setReportLevel(Log::error);
  EwaldRegion region(1, log);
  tools::Property prop = RegionDefinition(file);
  region.Initialize(prop.get("ewaldregion"));

  std::vector<PolarSegment> foreground = SingleProbe(5, probe_pos);
  region.ApplyFieldTo(foreground);

  const Eigen::Vector3d rhs = -foreground[0][0].V();
  const Eigen::Vector3d away = (probe_pos - charge_pos).normalized();
  const double alignment = rhs.normalized().dot(away);
  std::cout << "rhs (-V) = " << rhs.transpose()
            << "  alignment with 'away from +charge' = " << alignment
            << std::endl;

  BOOST_CHECK_GT(rhs.norm(), 0.0);
  // Comfortably positive rather than merely non-negative: a sign flip
  // gives -1 here, so anything near zero means the geometry of the test
  // is wrong rather than the code.
  BOOST_CHECK_GT(alignment, 0.9);

  std::remove(file.c_str());
}

// Flipping the background charge must flip the handed-over field exactly.
// The previous case pins the sign; this pins the linearity, so a
// magnitude error that survives a sign check is caught separately.
BOOST_AUTO_TEST_CASE(field_reverses_with_the_background_charge) {
  const Eigen::Vector3d charge_pos(10.0, 0.0, 0.0);
  const Eigen::Vector3d probe_pos = Eigen::Vector3d::Zero();

  auto field_for = [&](double q) {
    const std::string file = "ewaldregion_test_flip.hdf5";
    WriteSingleCharge(file, q, charge_pos, 5, probe_pos);
    Logger log;
    log.setReportLevel(Log::error);
    EwaldRegion region(1, log);
    tools::Property prop = RegionDefinition(file);
    region.Initialize(prop.get("ewaldregion"));
    std::vector<PolarSegment> foreground = SingleProbe(5, probe_pos);
    region.ApplyFieldTo(foreground);
    Eigen::Vector3d v = foreground[0][0].V();
    std::remove(file.c_str());
    return v;
  };

  const Eigen::Vector3d plus = field_for(+1.0);
  const Eigen::Vector3d minus = field_for(-1.0);
  BOOST_REQUIRE_GT(plus.norm(), 0.0);
  BOOST_CHECK_SMALL((plus + minus).norm() / plus.norm(), 1e-12);
}

// A foreground segment must interact with its OWN periodic images: it is
// carved out at ONE lattice site, and every other image of it is an
// ordinary background molecule. Legacy agrees (SetupMidground excludes
// by (segment id, na, nb, nc)); this code used to drop the segment at
// every translation, so those interactions vanished entirely.
//
// Production methane hid that completely -- neutral, and tetrahedral
// symmetry kills the dipole and traceless quadrupole, leaving an
// octupole that cancels by parity over a cubic image lattice. So the
// segment here is built with what methane lacks: neutral and dipole-free
// (the surface term then vanishes identically, since every piece of the
// shape bracket carries a Q0 or Q1 factor), but strongly quadrupolar.
//
// The reference is a direct lattice sum from Coulomb's law, absolutely
// convergent at 1/r^5. Before the fix the quantity is EXACTLY zero, so
// the discrimination does not rest on the tolerance.
BOOST_AUTO_TEST_CASE(foreground_segment_sees_its_own_periodic_images) {
  const std::string file = "ewaldregion_test_own_images.hdf5";
  const double box_length = 20.0;
  const double q = 0.5;
  const double d = 4.0;

  struct OwnImageSite {
    double charge;
    Eigen::Vector3d pos;
    Eigen::Vector3d dipole;
  };

  // Three sites along x. The charges (+q, -2q, +q) sum to zero with no
  // dipole by symmetry; the permanent dipoles sum to zero as well. The
  // cell carries neither a net charge nor a net dipole, so the direct
  // lattice sum below converges absolutely -- with a net dipole it is
  // conditionally convergent and its value depends on the shape of the
  // summation shell, a different question from the one asked here.
  const Eigen::Vector3d p0(-d, 0.0, 0.0);
  const Eigen::Vector3d p1 = Eigen::Vector3d::Zero();
  const Eigen::Vector3d p2(d, 0.0, 0.0);
  const Eigen::Vector3d m0(0.20, 0.00, 0.00);
  const Eigen::Vector3d m1(-0.10, 0.15, 0.00);
  const Eigen::Vector3d m2(-0.10, -0.15, 0.00);
  const Eigen::Vector3d no_dipole = Eigen::Vector3d::Zero();

  // RUN IN THREE CONFIGURATIONS. Charges alone was the original case,
  // and it cannot see anything wrong in the dipole path: with rank-0
  // sites getStaticDipole() is zero everywhere, so the mu term in the
  // structure factor, the dipole coincidence limits and dipolar Thole
  // are all multiplied away. Dipoles alone is the discriminating one --
  // any alpha-dependence it shows has to come from those branches, with
  // no charge channel to hide behind. The magnitudes are those of a
  // real rank-1 .mps (0.1-0.2 e*bohr), not the 5e-3 used elsewhere in
  // these tests, where a mu^2-sized error sits far below tolerance.
  const std::vector<std::pair<std::string, std::vector<OwnImageSite>>>
      configurations = {
          {"charges only",
           {{q, p0, no_dipole},
            {-2.0 * q, p1, no_dipole},
            {q, p2, no_dipole}}},
          {"dipoles only", {{0.0, p0, m0}, {0.0, p1, m1}, {0.0, p2, m2}}},
          {"charges and dipoles",
           {{q, p0, m0}, {-2.0 * q, p1, m1}, {q, p2, m2}}}};

  // SCANNED OVER ALPHA, and not for decoration. For the charge-only
  // configuration the erfc share of the own-image interaction, against
  // its 1.6e-3 hrt total, runs from 8.97e-19 at alpha = 0.50 to 1.03e-3
  // at alpha = 0.12. At 0.5 the reciprocal sum carries all of it, so a
  // single-alpha version tests only half the fix -- an earlier revision
  // did exactly that and passed with the real-space half reverted.
  auto ewald_energy_at = [&](const std::vector<OwnImageSite>& cluster,
                             double alpha) {
    {
      EwaldRegistry registry;
      PolarSegment seg("quad", 0);
      for (std::size_t i = 0; i < cluster.size(); ++i) {
        PolarSite site(Index(i), "C", cluster[i].pos);
        site.setpolarization(1e-6 * Eigen::Matrix3d::Identity());
        site.setCharge(cluster[i].charge);
        site.setStaticDipole(cluster[i].dipole);
        // No induced dipoles anywhere, so the returned energy is the
        // permanent-permanent channel alone.
        site.setInduced_Dipole(Eigen::Vector3d::Zero());
        seg.push_back(site);
      }
      registry.Register(0, EwaldChargeState::Neutral, seg);

      EwaldParameters params;
      params.alpha = alpha;
      // Held at a fixed multiple of alpha so each sum stays equally well
      // converged and the scan measures the splitting, not truncation.
      params.k_max = 12.0 * alpha;
      params.r_min = 12.0;
      params.field_tol = 1e-12;
      params.thole_a = 0.39;
      params.screening_factor = 6.0;
      params.shape = EwaldShape::Cube;
      params.box = box_length * Eigen::Matrix3d::Identity();

      CheckpointFile cpf(file, CheckpointAccessLevel::CREATE);
      CheckpointWriter w = cpf.getWriter();
      registry.WriteToCpt(w);
      CheckpointWriter wp = w.openChild("ewald_parameters");
      params.WriteToCpt(wp);
    }

    // The foreground: the same segment, at the same positions, carved out.
    PolarSegment fg("quad", 0);
    for (std::size_t i = 0; i < cluster.size(); ++i) {
      PolarSite site(Index(i), "C", cluster[i].pos);
      site.setpolarization(1e-6 * Eigen::Matrix3d::Identity());
      site.setCharge(cluster[i].charge);
      site.setStaticDipole(cluster[i].dipole);
      fg.push_back(site);
    }
    std::vector<PolarSegment> foreground{fg};

    Logger log;
    log.setReportLevel(Log::error);
    EwaldRegion region(1, log);
    tools::Property prop = RegionDefinition(file);
    region.Initialize(prop.get("ewaldregion"));
    const double e = region.ApplyFieldTo(foreground);
    std::remove(file.c_str());
    return e;
  };

  // Direct lattice sum: the cluster against every image of itself,
  // charge-charge plus charge-dipole plus dipole-dipole. The
  // charge-dipole term is written once per ordered pair and already
  // covers both orientations, so it is not halved.
  auto lattice_sum = [&](const std::vector<OwnImageSite>& cluster,
                         Index n_max) {
    double total = 0.0;
    for (Index na = -n_max; na <= n_max; ++na) {
      for (Index nb = -n_max; nb <= n_max; ++nb) {
        for (Index nc = -n_max; nc <= n_max; ++nc) {
          if (na == 0 && nb == 0 && nc == 0) {
            continue;
          }
          const Eigen::Vector3d L =
              box_length * Eigen::Vector3d(double(na), double(nb), double(nc));
          for (const auto& a : cluster) {
            for (const auto& b : cluster) {
              const Eigen::Vector3d r = a.pos - (b.pos + L);
              const double rn = r.norm();
              const double r3 = rn * rn * rn;
              total += a.charge * b.charge / rn;
              total += (a.charge * b.dipole.dot(r) -
                        b.charge * a.dipole.dot(r)) /
                       r3;
              total += a.dipole.dot(b.dipole) / r3 -
                       3.0 * a.dipole.dot(r) * b.dipole.dot(r) / (r3 * rn * rn);
            }
          }
        }
      }
    }
    return total;
  };

  for (const auto& configuration : configurations) {
    const std::vector<OwnImageSite>& cluster = configuration.second;

    const double reference_20 = lattice_sum(cluster, 20);
    const double reference_30 = lattice_sum(cluster, 30);
    const double tail = std::abs(reference_30 - reference_20);

    // Not vacuous: the interaction has to be a real number, far from the
    // zero the pre-fix code returns.
    BOOST_REQUIRE_GT(std::abs(reference_30), 1e-6);

    // The direct sum's own truncation is the limiting error here, not
    // the Ewald result, so the tolerance is taken from it rather than
    // guessed.
    const double tolerance =
        std::max(20.0 * tail, 1e-3 * std::abs(reference_30));

    std::cout << "own-image energy [" << configuration.first
              << "]: direct lattice sum = " << reference_30
              << " hrt, truncation tail = " << tail
              << ", tolerance = " << tolerance << std::endl;

    for (double alpha : {0.12, 0.20, 0.50}) {
      const double energy = ewald_energy_at(cluster, alpha);
      std::cout << "  alpha = " << alpha << "  Ewald = " << energy << " hrt"
                << std::endl;
      BOOST_CHECK_SMALL(std::abs(energy - reference_30), tolerance);
      // Stated separately: whatever the tolerance, the answer must not
      // be the zero the old convention returns.
      BOOST_CHECK_GT(std::abs(energy), 0.5 * std::abs(reference_30));
    }
  }
}

// PotentialAt against a direct lattice sum.
//
// A charge-only cluster alone in the cell: phi at each of its sites is
// the potential of every periodic image of itself, its own coincident
// copy suppressed. So sum_i q_i phi(r_i) IS the cluster's interaction
// with its own images -- the same quantity the direct sum above
// computes, and an independent reference rather than another route
// through the same code.
//
// Charges only on purpose: with dipoles the energy would need the field
// as well, and this case is about phi.
BOOST_AUTO_TEST_CASE(potential_at_reproduces_the_own_image_lattice_sum) {
  const std::string file = "ewaldregion_test_potential.hdf5";
  const double box_length = 20.0;
  const double q = 0.5;
  const double d = 4.0;

  const std::vector<std::pair<double, Eigen::Vector3d>> cluster = {
      {q, Eigen::Vector3d(-d, 0.0, 0.0)},
      {-2.0 * q, Eigen::Vector3d::Zero()},
      {q, Eigen::Vector3d(d, 0.0, 0.0)}};

  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  for (const auto& entry : cluster) {
    centroid += entry.second;
  }
  centroid /= double(cluster.size());

  auto potential_energy_at = [&](double alpha) {
    {
      EwaldRegistry registry;
      PolarSegment seg("quad", 0);
      for (std::size_t i = 0; i < cluster.size(); ++i) {
        PolarSite site(Index(i), "C", cluster[i].second);
        site.setpolarization(1e-6 * Eigen::Matrix3d::Identity());
        site.setCharge(cluster[i].first);
        site.setInduced_Dipole(Eigen::Vector3d::Zero());
        seg.push_back(site);
      }
      registry.Register(0, EwaldChargeState::Neutral, seg);

      EwaldParameters params;
      params.alpha = alpha;
      params.k_max = 12.0 * alpha;
      params.r_min = 12.0;
      params.field_tol = 1e-12;
      params.thole_a = 0.39;
      params.screening_factor = 6.0;
      params.shape = EwaldShape::Cube;
      params.box = box_length * Eigen::Matrix3d::Identity();

      CheckpointFile cpf(file, CheckpointAccessLevel::CREATE);
      CheckpointWriter w = cpf.getWriter();
      registry.WriteToCpt(w);
      CheckpointWriter wp = w.openChild("ewald_parameters");
      params.WriteToCpt(wp);
    }

    Logger log;
    log.setReportLevel(Log::error);
    EwaldRegion region(1, log);
    tools::Property prop = RegionDefinition(file);
    region.Initialize(prop.get("ewaldregion"));

    // A point carries no segment identity, so the foreground has to be
    // declared rather than inferred from a call argument.
    region.RegisterForeground({{0, centroid}});

    std::vector<Eigen::Vector3d> points;
    for (const auto& entry : cluster) {
      points.push_back(entry.second);
    }
    const Eigen::VectorXd phi = region.PotentialAt(points);
    std::remove(file.c_str());

    BOOST_REQUIRE_EQUAL(phi.size(), Index(cluster.size()));
    double energy = 0.0;
    for (std::size_t i = 0; i < cluster.size(); ++i) {
      energy += cluster[i].first * phi[Index(i)];
    }
    return energy;
  };

  auto lattice_sum = [&](Index n_max) {
    double total = 0.0;
    for (Index na = -n_max; na <= n_max; ++na) {
      for (Index nb = -n_max; nb <= n_max; ++nb) {
        for (Index nc = -n_max; nc <= n_max; ++nc) {
          if (na == 0 && nb == 0 && nc == 0) {
            continue;
          }
          const Eigen::Vector3d L =
              box_length * Eigen::Vector3d(double(na), double(nb), double(nc));
          for (const auto& a : cluster) {
            for (const auto& b : cluster) {
              total += a.first * b.first / (a.second - (b.second + L)).norm();
            }
          }
        }
      }
    }
    return total;
  };

  const double reference_20 = lattice_sum(20);
  const double reference_30 = lattice_sum(30);
  const double tail = std::abs(reference_30 - reference_20);
  const double tolerance =
      std::max(20.0 * tail, 1e-3 * std::abs(reference_30));

  BOOST_REQUIRE_GT(std::abs(reference_30), 1e-6);
  std::cout << "PotentialAt own-image energy: direct lattice sum = "
            << reference_30 << " hrt, tolerance = " << tolerance << std::endl;

  // Scanned, for the same reason the field case is: alpha moves weight
  // between the real and reciprocal halves of phi, and only the total is
  // meaningful.
  for (double alpha : {0.12, 0.20, 0.50}) {
    const double energy = potential_energy_at(alpha);
    std::cout << "  alpha = " << alpha << "  from PotentialAt = " << energy
              << " hrt" << std::endl;
    BOOST_CHECK_SMALL(std::abs(energy - reference_30), tolerance);
    BOOST_CHECK_GT(std::abs(energy), 0.5 * std::abs(reference_30));
  }
}

// PotentialAt and ApplyFieldTo, held against each other.
//
// These are the two routes a job takes to the same background: a QM
// region asks for phi on a grid, a polar region asks for the field and
// energy on its sites. They share a registry and a foreground
// declaration and nothing else -- PotentialAt is cache-free, duplicates
// the distance cull, and does its own erf bookkeeping, precisely because
// a DFT grid cannot be walked with the neighbour cache the site path
// relies on. Duplicated rules drift, and until now nothing compared the
// two on real multipoles.
//
// For a RANK-0 foreground the comparison is exact rather than
// approximate. ApplyFieldTo's returned energy is then sum_i q_i phi(r_i)
// and nothing else, since every dipole-field term it also assembles
// contracts against a permanent dipole that is zero. PotentialAt returns
// that same phi at the same points. So one contracted against the
// foreground's charges must reproduce the other to round-off.
//
// NO INDUCED DIPOLES in the background here, deliberately. The induced
// channel is the one place these two legitimately differ: ApplyFieldTo's
// target is a polarizable site and is Thole-damped against the
// background's induced dipoles, while PotentialAt's target is a point in
// space, which has no polarizability to damp with (see
// probe_is_not_thole_damped_against_an_induced_dipole in
// test_ewaldrealspacesum). Letting that in would make any disagreement
// here ambiguous. This case is about the permanent channel, where the
// two have no freedom to differ at all.
BOOST_AUTO_TEST_CASE(potential_at_matches_apply_field_to_on_a_rank0_foreground) {
  const std::string file = "ewaldregion_test_crosscheck.hdf5";
  const double box_length = 24.0;
  const double alpha = 0.35;

  // Three neutral rank-0 segments, one of which becomes the foreground.
  // Neutral per segment so the k = 0 omission carries no net charge.
  const std::vector<std::vector<std::pair<double, Eigen::Vector3d>>> segments = {
      {{-0.6, Eigen::Vector3d(0.0, 0.0, 0.0)},
       {0.3, Eigen::Vector3d(1.4, 0.2, -0.3)},
       {0.3, Eigen::Vector3d(-1.1, 0.9, 0.5)}},
      {{-0.5, Eigen::Vector3d(7.0, 1.0, 2.0)},
       {0.5, Eigen::Vector3d(8.3, 1.6, 2.4)}},
      {{-0.4, Eigen::Vector3d(3.0, 6.5, -4.0)},
       {0.4, Eigen::Vector3d(4.1, 7.2, -3.4)}}};

  auto centroid_of = [](const std::vector<std::pair<double, Eigen::Vector3d>>&
                            sites) {
    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    for (const auto& s : sites) {
      c += s.second;
    }
    return Eigen::Vector3d(c / double(sites.size()));
  };

  {
    EwaldRegistry registry;
    for (std::size_t s = 0; s < segments.size(); ++s) {
      PolarSegment seg("seg", Index(s));
      for (std::size_t i = 0; i < segments[s].size(); ++i) {
        PolarSite site(Index(i), "C", segments[s][i].second);
        site.setpolarization(1e-6 * Eigen::Matrix3d::Identity());
        site.setCharge(segments[s][i].first);
        // Explicitly zero: this case compares the permanent channel.
        site.setInduced_Dipole(Eigen::Vector3d::Zero());
        seg.push_back(site);
      }
      registry.Register(Index(s), EwaldChargeState::Neutral, seg);
    }

    EwaldParameters params;
    params.alpha = alpha;
    params.k_max = 12.0 * alpha;
    params.r_min = 12.0;
    params.field_tol = 1e-12;
    params.thole_a = 0.39;
    params.screening_factor = 6.0;
    params.shape = EwaldShape::Cube;
    params.box = box_length * Eigen::Matrix3d::Identity();

    CheckpointFile cpf(file, CheckpointAccessLevel::CREATE);
    CheckpointWriter w = cpf.getWriter();
    registry.WriteToCpt(w);
    CheckpointWriter wp = w.openChild("ewald_parameters");
    params.WriteToCpt(wp);
  }

  Logger log;
  log.setReportLevel(Log::error);
  EwaldRegion region(1, log);
  tools::Property prop = RegionDefinition(file);
  region.Initialize(prop.get("ewaldregion"));

  const Index fg_id = 0;
  region.RegisterForeground({{fg_id, centroid_of(segments[fg_id])}});

  // The foreground as a polar region would hand it over: the same sites,
  // at the same positions, with the same charges.
  std::vector<PolarSegment> foreground;
  {
    PolarSegment seg("seg", fg_id);
    for (std::size_t i = 0; i < segments[fg_id].size(); ++i) {
      PolarSite site(Index(i), "C", segments[fg_id][i].second);
      site.setpolarization(1e-6 * Eigen::Matrix3d::Identity());
      site.setCharge(segments[fg_id][i].first);
      site.setInduced_Dipole(Eigen::Vector3d::Zero());
      seg.push_back(site);
    }
    foreground.push_back(seg);
  }

  // ApplyFieldTo first: it builds the sums, and PotentialAt then reuses
  // the very same objects. Sharing them is the point -- a disagreement
  // here can only come from the two traversals, not from two different
  // backgrounds.
  const double energy_from_field_path = region.ApplyFieldTo(foreground);

  std::vector<Eigen::Vector3d> points;
  for (const auto& s : segments[fg_id]) {
    points.push_back(s.second);
  }
  const Eigen::VectorXd phi = region.PotentialAt(points);
  std::remove(file.c_str());

  BOOST_REQUIRE_EQUAL(phi.size(), Index(points.size()));
  double energy_from_potential = 0.0;
  for (std::size_t i = 0; i < segments[fg_id].size(); ++i) {
    energy_from_potential += segments[fg_id][i].first * phi[Index(i)];
  }

  // Not vacuous: a background this close produces a real interaction, so
  // agreeing on zero would not be agreeing on anything.
  BOOST_REQUIRE_GT(std::abs(energy_from_field_path), 1e-6);

  BOOST_CHECK_CLOSE(energy_from_potential, energy_from_field_path, 1e-8);
}

// The same cross-check as above, with a CHARGED foreground -- which is
// what tests the gauge, and what the case above cannot.
//
// phi carries an arbitrary additive constant: the k = 0 term is omitted
// (the uniform neutralising background), so phi is only defined up to
// phi_0. A region of net charge q embedded in it shifts by q * phi_0.
// The neutral case above therefore says nothing about phi_0 at all --
// sum q_i = 0 makes it drop out identically, however wrong it might be.
//
// That constant does NOT cancel in a site energy, which is the whole
// point of this code. The neutral job has q_QM = 0 and the hole job has
// q_QM = +1, so E(h) - E(n) carries +1 * phi_0 outright.
//
// So: background stays neutral (the registry is unchanged, and the
// suppressed copy is still the neutral one the reciprocal sum actually
// placed), and only the foreground handed to ApplyFieldTo carries net
// charge -- exactly the arrangement a hole job produces, where the
// background holds a segment's neutral copy and the job replaces it with
// a charged one. sum_i q_i phi(r_i) then weights phi_0 by sum q_i, and a
// disagreement between the two paths proportional to that sum is the one
// thing nothing else in this file can produce.
//
// What this does and does not establish: it shows the two paths share a
// gauge, not that the shared gauge is the right one. The chain closes
// elsewhere -- ApplyFieldTo's convention is the one validated against
// legacy by the MM/MM runs, so PotentialAt agreeing with it here
// inherits that validation. Which is what a QM/MM site energy needs,
// since its QM half arrives through PotentialAt and its polar half
// through ApplyFieldTo, in the same job.
BOOST_AUTO_TEST_CASE(potential_at_matches_apply_field_to_on_a_charged_foreground) {
  const std::string file = "ewaldregion_test_crosscheck_charged.hdf5";
  const double box_length = 24.0;
  const double alpha = 0.35;

  // Neutral background, as a converged background always is.
  const std::vector<std::vector<std::pair<double, Eigen::Vector3d>>> background =
      {{{-0.6, Eigen::Vector3d(0.0, 0.0, 0.0)},
        {0.3, Eigen::Vector3d(1.4, 0.2, -0.3)},
        {0.3, Eigen::Vector3d(-1.1, 0.9, 0.5)}},
       {{-0.5, Eigen::Vector3d(7.0, 1.0, 2.0)},
        {0.5, Eigen::Vector3d(8.3, 1.6, 2.4)}},
       {{-0.4, Eigen::Vector3d(3.0, 6.5, -4.0)},
        {0.4, Eigen::Vector3d(4.1, 7.2, -3.4)}}};

  // The job's charge state for segment 0: same sites, same positions,
  // charges summing to +1 rather than 0. Deliberately not a uniform
  // shift of the neutral values, so no accidental symmetry can make the
  // comparison pass for the wrong reason.
  const std::vector<double> charged_q = {0.5, 0.2, 0.3};

  const Index fg_id = 0;
  BOOST_REQUIRE_EQUAL(charged_q.size(), background[fg_id].size());

  double net_charge = 0.0;
  for (double q : charged_q) {
    net_charge += q;
  }
  // The case is about phi_0's weight. If this were zero it would silently
  // become the neutral case again and test nothing new.
  BOOST_REQUIRE_GT(std::abs(net_charge), 0.5);

  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  for (const auto& s : background[fg_id]) {
    centroid += s.second;
  }
  centroid /= double(background[fg_id].size());

  {
    EwaldRegistry registry;
    for (std::size_t s = 0; s < background.size(); ++s) {
      PolarSegment seg("seg", Index(s));
      for (std::size_t i = 0; i < background[s].size(); ++i) {
        PolarSite site(Index(i), "C", background[s][i].second);
        site.setpolarization(1e-6 * Eigen::Matrix3d::Identity());
        site.setCharge(background[s][i].first);
        site.setInduced_Dipole(Eigen::Vector3d::Zero());
        seg.push_back(site);
      }
      registry.Register(Index(s), EwaldChargeState::Neutral, seg);
    }

    EwaldParameters params;
    params.alpha = alpha;
    params.k_max = 12.0 * alpha;
    params.r_min = 12.0;
    params.field_tol = 1e-12;
    params.thole_a = 0.39;
    params.screening_factor = 6.0;
    params.shape = EwaldShape::Cube;
    params.box = box_length * Eigen::Matrix3d::Identity();

    CheckpointFile cpf(file, CheckpointAccessLevel::CREATE);
    CheckpointWriter w = cpf.getWriter();
    registry.WriteToCpt(w);
    CheckpointWriter wp = w.openChild("ewald_parameters");
    params.WriteToCpt(wp);
  }

  Logger log;
  log.setReportLevel(Log::error);
  EwaldRegion region(1, log);
  tools::Property prop = RegionDefinition(file);
  region.Initialize(prop.get("ewaldregion"));
  region.RegisterForeground({{fg_id, centroid}});

  std::vector<PolarSegment> foreground;
  {
    PolarSegment seg("seg", fg_id);
    for (std::size_t i = 0; i < background[fg_id].size(); ++i) {
      PolarSite site(Index(i), "C", background[fg_id][i].second);
      site.setpolarization(1e-6 * Eigen::Matrix3d::Identity());
      site.setCharge(charged_q[i]);
      site.setInduced_Dipole(Eigen::Vector3d::Zero());
      seg.push_back(site);
    }
    foreground.push_back(seg);
  }

  // ApplyFieldTo first, so PotentialAt reuses the very same sums.
  const double energy_from_field_path = region.ApplyFieldTo(foreground);

  std::vector<Eigen::Vector3d> points;
  for (const auto& s : background[fg_id]) {
    points.push_back(s.second);
  }
  const Eigen::VectorXd phi = region.PotentialAt(points);
  std::remove(file.c_str());

  BOOST_REQUIRE_EQUAL(phi.size(), Index(points.size()));
  double energy_from_potential = 0.0;
  for (std::size_t i = 0; i < charged_q.size(); ++i) {
    energy_from_potential += charged_q[i] * phi[Index(i)];
  }

  BOOST_REQUIRE_GT(std::abs(energy_from_field_path), 1e-6);
  BOOST_CHECK_CLOSE(energy_from_potential, energy_from_field_path, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
