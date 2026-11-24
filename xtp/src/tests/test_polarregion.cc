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

#define BOOST_TEST_MODULE polarregion_test

// Standard includes
#include <memory>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/polarregion.h"
#include "votca/xtp/staticregion.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(polarregion_test)

namespace {

// A small polar molecule: a charged pair with a permanent dipole on one end,
// isotropically polarizable, centred at `c` (bohr).
PolarSegment Molecule(Index id, const Eigen::Vector3d& c) {
  PolarSegment seg("mol", id);
  PolarSite a(0, "C", c + Eigen::Vector3d(0.0, 0.0, -1.1));
  a.setCharge(0.3);
  a.setStaticDipole(Eigen::Vector3d(0.05, -0.02, 0.1));
  a.setpolarization(8.0 * Eigen::Matrix3d::Identity());
  PolarSite b(1, "O", c + Eigen::Vector3d(0.0, 0.0, 1.1));
  b.setCharge(-0.3);
  b.setpolarization(5.0 * Eigen::Matrix3d::Identity());
  seg.push_back(a);
  seg.push_back(b);
  return seg;
}

std::vector<PolarSegment> Molecules() {
  std::vector<PolarSegment> mols;
  Index id = 0;
  for (Index i = 0; i < 3; ++i) {
    for (Index j = 0; j < 2; ++j) {
      const Eigen::Vector3d c(7.0 * double(i), 7.5 * double(j),
                              0.6 * double(i * j));
      mols.push_back(Molecule(id++, c));
    }
  }
  return mols;
}

// Fixed external charges: the field the polar regions respond to.
std::unique_ptr<StaticRegion> External(Index id, Logger& log) {
  auto region = std::make_unique<StaticRegion>(id, log);
  StaticSegment seg("ext", 100);
  StaticSite p(0, "Na", Eigen::Vector3d(7.0, 3.7, 9.0));
  p.setCharge(1.0);
  StaticSite m(1, "Cl", Eigen::Vector3d(-6.0, 3.7, -2.0));
  m.setCharge(-0.7);
  seg.push_back(p);
  seg.push_back(m);
  region->push_back(seg);
  return region;
}

std::unique_ptr<PolarRegion> Polar(Index id, Logger& log,
                                   const std::vector<PolarSegment>& segs,
                                   const std::string& guess = "") {
  auto region = std::make_unique<PolarRegion>(id, log);
  for (const PolarSegment& s : segs) {
    region->push_back(s);
  }
  votca::tools::Property opt;
  opt.add("max_iter", "200");
  opt.add("tolerance_dipole", "1e-12");
  opt.add("tolerance_energy", "1e-12");
  opt.add("exp_damp", "0.39");
  if (!guess.empty()) {
    opt.add("guess", guess);
  }
  region->Initialize(opt);
  return region;
}

// The inter-region loop of the QMMM calculator, run to convergence.
double Converge(std::vector<std::unique_ptr<Region>>& regions) {
  double etot = 0.0;
  for (Index iter = 0; iter < 60; ++iter) {
    for (std::unique_ptr<Region>& r : regions) {
      r->Reset();
      r->Evaluate(regions);
      BOOST_REQUIRE(r->Successful());
    }
    etot = 0.0;
    for (const std::unique_ptr<Region>& r : regions) {
      etot += r->Etotal();
    }
  }
  return etot;
}

std::vector<Eigen::Vector3d> Dipoles(
    const std::vector<std::unique_ptr<Region>>& regions) {
  std::vector<Eigen::Vector3d> out;
  for (const std::unique_ptr<Region>& r : regions) {
    const auto* polar = dynamic_cast<const PolarRegion*>(r.get());
    if (polar == nullptr) {
      continue;
    }
    for (const PolarSegment& seg : *polar) {
      for (const PolarSite& site : seg) {
        out.push_back(site.Induced_Dipole());
      }
    }
  }
  return out;
}

}  // namespace

// The same polar molecules in the same external field, once as one polar
// region and once split into two, must give the same induced dipoles and
// the same total energy: the partition into regions is bookkeeping only.
// Before the fix the inter-region static energy of the lower-id region was
// summed triangularly over the other region's segments, which made the
// split total wrong and, as it depends on the other region's induced
// dipoles, not stationary -- visible as a 0.5% error in QM/MM forces with
// two polarizable regions.
BOOST_AUTO_TEST_CASE(split_equals_merged) {
  Logger log;
  log.setReportLevel(votca::Log::error);
  const std::vector<PolarSegment> mols = Molecules();

  std::vector<std::unique_ptr<Region>> merged;
  merged.push_back(Polar(0, log, mols));
  merged.push_back(External(1, log));
  const double e_merged = Converge(merged);

  const std::vector<PolarSegment> inner(mols.begin(), mols.begin() + 2);
  const std::vector<PolarSegment> outer(mols.begin() + 2, mols.end());
  std::vector<std::unique_ptr<Region>> split;
  split.push_back(Polar(0, log, inner));
  split.push_back(Polar(1, log, outer));
  split.push_back(External(2, log));
  const double e_split = Converge(split);

  const std::vector<Eigen::Vector3d> mu_merged = Dipoles(merged);
  const std::vector<Eigen::Vector3d> mu_split = Dipoles(split);
  BOOST_REQUIRE_EQUAL(mu_merged.size(), mu_split.size());
  double dmu = 0.0;
  double mu = 0.0;
  for (std::size_t i = 0; i < mu_merged.size(); ++i) {
    dmu = std::max(dmu, (mu_merged[i] - mu_split[i]).norm());
    mu = std::max(mu, mu_merged[i].norm());
  }
  std::cout << "E merged " << e_merged << " E split " << e_split << " diff "
            << e_split - e_merged << "  max|dmu| " << dmu << " of " << mu
            << std::endl;
  BOOST_CHECK_GT(mu, 1e-3);
  BOOST_CHECK_LT(dmu, 1e-8);
  BOOST_CHECK_CLOSE(e_split, e_merged, 1e-6);
}

// polarregion guess: the induced moments of a checkpoint of an earlier job
// (here: the converged merged region) start the first solve. Segments are
// matched by id, since the checkpoint does not keep their order; a
// checkpoint whose segments do not match is ignored.
BOOST_AUTO_TEST_CASE(guess_from_checkpoint) {
  Logger log;
  log.setReportLevel(votca::Log::error);
  const std::vector<PolarSegment> mols = Molecules();
  std::vector<std::unique_ptr<Region>> converged;
  converged.push_back(Polar(0, log, mols));
  converged.push_back(External(1, log));
  Converge(converged);
  const std::vector<Eigen::Vector3d> mu = Dipoles(converged);
  const std::string file = "polarregion_guess.hdf5";
  {
    CheckpointFile cpf(file, CheckpointAccessLevel::CREATE);
    CheckpointWriter w = cpf.getWriter("region_0");
    converged[0]->WriteToCpt(w);
  }
  // the same segments: moments loaded before any evaluation
  std::vector<std::unique_ptr<Region>> seeded;
  seeded.push_back(Polar(0, log, mols, file));
  const std::vector<Eigen::Vector3d> mu_seed = Dipoles(seeded);
  BOOST_REQUIRE_EQUAL(mu_seed.size(), mu.size());
  double dev = 0.0;
  for (std::size_t i = 0; i < mu.size(); ++i) {
    dev = std::max(dev, (mu_seed[i] - mu[i]).norm());
  }
  BOOST_CHECK_SMALL(dev, 1e-12);
  // and converging from them gives the same answer
  seeded.push_back(External(1, log));
  const double e_seeded = Converge(seeded);
  const double e_ref = converged[0]->Etotal() + converged[1]->Etotal();
  BOOST_CHECK_CLOSE(e_seeded, e_ref, 1e-8);
  // different segments (one molecule fewer): ignored, moments stay zero
  const std::vector<PolarSegment> fewer(mols.begin(), mols.end() - 1);
  std::vector<std::unique_ptr<Region>> other;
  other.push_back(Polar(0, log, fewer, file));
  for (const Eigen::Vector3d& d : Dipoles(other)) {
    BOOST_CHECK_EQUAL(d.norm(), 0.0);
  }
  std::remove(file.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
