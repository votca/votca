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
#include <fstream>
#include <iomanip>
#include <numeric>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/statesaver.h"

// Local private VOTCA includes
#include "latticegen.h"

namespace votca {
namespace xtp {

void LatticeGen::ParseOptions(const tools::Property& options) {
  statefile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".statefile", job_name_ + ".hdf5");
  energyfile_ =
      options.ifExistsReturnElseReturnDefault<std::string>(".energyfile", "");

  std::vector<Index> dims = options.get(".dims").as<std::vector<Index>>();
  if (dims.size() != 3) {
    throw std::runtime_error("latticegen: dims needs three numbers");
  }
  lattice_.dims = {dims[0], dims[1], dims[2]};
  lattice_.spacing =
      options.get(".spacing").as<double>() * tools::conv::nm2bohr;
  lattice_.cutoff = options.get(".cutoff").as<double>() * tools::conv::nm2bohr;
  double ev2 = tools::conv::ev2hrt * tools::conv::ev2hrt;
  lattice_.J2 = options.get(".J2").as<double>() * ev2;
  lattice_.decay_length =
      options.get(".decaylength").as<double>() * tools::conv::nm2bohr;
  lattice_.reorganization =
      options.get(".reorganization").as<double>() * tools::conv::ev2hrt;
  lattice_.mean_energy =
      options.get(".energy.mean").as<double>() * tools::conv::ev2hrt;
  lattice_.sigma =
      options.get(".energy.sigma").as<double>() * tools::conv::ev2hrt;
  lattice_.disorder = LatticeModel::DisorderFromString(
      options.get(".energy.disorder").as<std::string>());
  lattice_.eta_a = options.get(".energy.eta_a").as<double>();
  lattice_.eta_b = options.get(".energy.eta_b").as<double>();
  lattice_.l_b = options.get(".energy.l_b").as<Index>();
  lattice_.l_c = options.get(".energy.l_c").as<Index>();
  lattice_.include_self = options.get(".energy.includeself").as<bool>();
  lattice_.seed = options.get(".energy.seed").as<std::uint64_t>();

  std::vector<std::string> states =
      options.get(".states").as<std::vector<std::string>>();
  lattice_.states.clear();
  for (const std::string& s : states) {
    QMStateType type(s);
    if (!type.isKMCState()) {
      throw std::runtime_error("latticegen: " + s +
                               " is not a hopping state (e, h, s, t)");
    }
    lattice_.states.push_back(type);
  }
  lattice_.segment_name = options.get(".segmentname").as<std::string>();
}

bool LatticeGen::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  if (tools::filesystem::FileExists(statefile_)) {
    throw std::runtime_error("latticegen: " + statefile_ +
                             " exists already, remove it or choose another "
                             "statefile name");
  }

  LatticeModel model(lattice_);
  Topology top;
  model.BuildTopology(top);

  const std::vector<double>& e = model.SiteEnergies();
  double mean = std::accumulate(e.begin(), e.end(), 0.0) / double(e.size());
  double var = std::accumulate(e.begin(), e.end(), 0.0,
                               [mean](double acc, double x) {
                                 return acc + (x - mean) * (x - mean);
                               }) /
               double(e.size());
  XTP_LOG(Log::error, log_)
      << "Lattice " << lattice_.dims[0] << "x" << lattice_.dims[1] << "x"
      << lattice_.dims[2] << ": " << model.NumberOfSites() << " sites, "
      << top.NBList().size() << " pairs" << std::flush;
  XTP_LOG(Log::error, log_)
      << "Site energies: mean " << mean * tools::conv::hrt2ev
      << " eV, standard deviation " << std::sqrt(var) * tools::conv::hrt2ev
      << " eV" << std::flush;

  StateSaver saver(statefile_);
  saver.WriteFrame(top);
  XTP_LOG(Log::error, log_) << "Wrote " << statefile_ << std::flush;

  if (!energyfile_.empty()) {
    std::ofstream out(energyfile_);
    out << "# segment id, x, y, z [nm], site energy [eV]\n";
    out << std::setprecision(12);
    const std::vector<Eigen::Vector3d>& pos = model.Positions();
    for (Index i = 0; i < model.NumberOfSites(); ++i) {
      Eigen::Vector3d p = pos[i] * tools::conv::bohr2nm;
      out << i << " " << p.x() << " " << p.y() << " " << p.z() << " "
          << e[i] * tools::conv::hrt2ev << "\n";
    }
    XTP_LOG(Log::error, log_)
        << "Wrote site energies to " << energyfile_ << std::flush;
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
