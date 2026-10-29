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

// Standard includes
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/filesystem.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/statesaver.h"

// Local private VOTCA includes
#include "networkgen.h"

namespace votca {
namespace xtp {

void NetworkGen::ParseOptions(const tools::Property& options) {
  statefile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".statefile", job_name_ + ".hdf5");
  energyfile_ =
      options.ifExistsReturnElseReturnDefault<std::string>(".energyfile", "");
  parameterfile_ = options.get(".parameters").as<std::string>();
  sites_ = options.ifExistsReturnElseReturnDefault<Index>(".sites", 0);
  if (options.exists(".box") &&
      !options.get(".box").as<std::string>().empty()) {
    box_ = options.get(".box").as<std::vector<double>>();
    if (box_.size() != 3) {
      throw std::runtime_error("networkgen: box needs three lengths");
    }
  }
  seed_ = options.ifExistsReturnElseReturnDefault<Index>(".seed", 0);
  temperature_ = options.get(".temperature").as<double>();
  trap_factor_ = options.get(".trapfactor").as<double>();
}

bool NetworkGen::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  if (tools::filesystem::FileExists(statefile_)) {
    throw std::runtime_error("networkgen: " + statefile_ +
                             " exists already, remove it or choose another "
                             "statefile name");
  }
  tools::Property file;
  file.LoadFromXML(parameterfile_);
  StochasticNetwork::Options opt =
      StochasticNetwork::FromProperty(file.get("stochasticnetwork"));
  const Index fitted_sites = opt.sites;
  const Eigen::Vector3d fitted_box = opt.box;
  if (sites_ > 0) {
    opt.sites = sites_;
  }
  if (!box_.empty()) {
    opt.box = Eigen::Vector3d(box_[0], box_[1], box_[2]) * tools::conv::nm2bohr;
  } else if (fitted_sites > 0 && fitted_box.minCoeff() > 0.0) {
    // the density and aspect ratio of the fitted system
    opt.box = fitted_box * std::cbrt(double(opt.sites) / double(fitted_sites));
  } else {
    throw std::runtime_error(
        "networkgen: give box, the parameter file has no box and sites");
  }
  if (seed_ > 0) {
    opt.seed = std::uint64_t(seed_);
  }
  opt.trap_factor = trap_factor_;
  Eigen::Vector3d box_nm = opt.box * tools::conv::bohr2nm;
  XTP_LOG(Log::error, log_)
      << "Stochastic network from " << parameterfile_ << ": " << opt.sites
      << " sites in a box of " << box_nm.x() << " x " << box_nm.y() << " x "
      << box_nm.z() << " nm^3, seed " << opt.seed << std::flush;

  StochasticNetwork model(opt);
  Topology top;
  model.BuildTopology(top);

  XTP_LOG(Log::error, log_)
      << "Point process: " << model.Iterations() << " iterations; "
      << top.NBList().size() << " pairs, "
      << 2.0 * double(top.NBList().size()) / double(model.NumberOfSites())
      << " neighbours per site" << std::flush;
  if (!opt.traps.empty()) {
    std::ostringstream t;
    t << "Traps: " << model.Traps().size() << " inserted (" << opt.traps.size()
      << " template(s), densities x " << trap_factor_ << ")";
    for (std::size_t k = 0; k < opt.traps.size(); ++k) {
      const auto& tr = opt.traps[k];
      Index count =
          Index(std::count_if(model.Traps().begin(), model.Traps().end(),
                              [k](const StochasticNetwork::InsertedTrap& x) {
                                return x.type == Index(k);
                              }));
      t << "\n    template " << k << ": "
        << opt.carriers[tr.carrier].state.ToString() << " trap of "
        << tr.offsets.size() << " site(s), deepest "
        << tr.offsets[0][tr.carrier] * tools::conv::hrt2ev
        << " eV from the bulk mean, density " << tr.density * trap_factor_
        << " per site: " << count << " inserted";
    }
    XTP_LOG(Log::error, log_) << t.str() << std::flush;
  }
  const double kT = tools::conv::kB * temperature_;  // eV
  for (std::size_t c = 0; c < opt.carriers.size(); ++c) {
    const std::vector<double>& e = model.SiteEnergies(Index(c));
    double mean = std::accumulate(e.begin(), e.end(), 0.0) / double(e.size());
    double var = std::accumulate(e.begin(), e.end(), 0.0,
                                 [mean](double acc, double x) {
                                   return acc + (x - mean) * (x - mean);
                                 }) /
                 double(e.size());
    double sigma = std::sqrt(var) * tools::conv::hrt2ev;
    // Gaussian DOS: occupied states at mean - sigma^2/kT need
    // (sigma/kT)^2 <= -5.7 + 1.05 ln N (Baumeier et al., PRB 86, 184202)
    double s2 = (sigma / kT) * (sigma / kT);
    double nmin = std::exp((s2 + 5.7) / 1.05);
    XTP_LOG(Log::error, log_)
        << "Carrier " << opt.carriers[c].state.ToString()
        << ": site energies mean " << mean * tools::conv::hrt2ev
        << " eV, standard deviation " << sigma << " eV = " << sigma / kT
        << " kT at " << temperature_ << " K; the density of occupied states "
        << (double(opt.sites) >= nmin ? "is" : "is NOT")
        << " sampled (needs about " << std::setprecision(2) << std::scientific
        << nmin << " sites)" << std::defaultfloat << std::setprecision(6)
        << std::flush;
  }

  StateSaver saver(statefile_);
  saver.WriteFrame(top);
  XTP_LOG(Log::error, log_) << "Wrote " << statefile_ << std::flush;

  if (!energyfile_.empty()) {
    std::ofstream out(energyfile_);
    out << "# segment id, x, y, z [nm], site energies [eV] of";
    for (const auto& c : opt.carriers) {
      out << " " << c.state.ToString();
    }
    out << "\n" << std::setprecision(12);
    const std::vector<Eigen::Vector3d>& pos = model.Positions();
    for (Index i = 0; i < model.NumberOfSites(); ++i) {
      Eigen::Vector3d p = pos[i] * tools::conv::bohr2nm;
      out << i << " " << p.x() << " " << p.y() << " " << p.z();
      for (std::size_t c = 0; c < opt.carriers.size(); ++c) {
        out << " " << model.SiteEnergies(Index(c))[i] * tools::conv::hrt2ev;
      }
      out << "\n";
    }
    XTP_LOG(Log::error, log_)
        << "Wrote site energies to " << energyfile_ << std::flush;
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
