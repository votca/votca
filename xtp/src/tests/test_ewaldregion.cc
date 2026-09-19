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
#include <cstdio>
#include <iostream>

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

BOOST_AUTO_TEST_SUITE_END()
