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

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/resultfile.h"
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "grw.h"

namespace votca {
namespace xtp {

void GRW::ParseOptions(const tools::Property& options) {
  carriertype_ = QMStateType(options.get(".carriertype").as<std::string>());
  if (!carriertype_.isKMCState()) {
    throw std::runtime_error("grw cannot be run for state: " +
                             carriertype_.ToLongString());
  }
  temperature_ = options.get(".temperature").as<double>() * tools::conv::kB *
                 tools::conv::ev2hrt;
  field_Vm_ = options.get(".field").as<Eigen::Vector3d>();
  // V/m -> Hartree/bohr
  field_ = field_Vm_ * (tools::conv::ev2hrt / (1e9 * tools::conv::nm2bohr));
  tof_ = ToFGraph::OptionsFromProperty(options.get(".tof"));
  carriers_ = options.get(".numberofcarriers").as<Index>();
  observable_ = GraphRandomWalk::ObservableFromString(
      options.get(".observable").as<std::string>());
  settings_.max_states = options.get(".maxstates").as<Index>();
  solver_ = GraphRandomWalk::SolverFromString(
      options.get(".solver").as<std::string>());
  settings_.direct_limit = options.get(".directlimit").as<Index>();
  settings_.tolerance = options.get(".tolerance").as<double>();
  settings_.max_iterations = options.get(".maxiterations").as<Index>();
  sourcefile_ = options.get(".sourcefile").as<std::string>();
  occupationfile_ = options.get(".occupationfile").as<std::string>();
  resultfile_ = options.get(".resultfile").as<std::string>();
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n ...");
}

bool GRW::Evaluate(Topology& top) {
  XTP_LOG(Log::error, log_) << "\n-----------------------------------"
                               "\n  GRAPH RANDOM WALK TIME OF FLIGHT"
                               "\n-----------------------------------"
                            << std::flush;
  if (carriers_ != 1 && carriers_ != 2) {
    throw std::runtime_error(
        "grw: numberofcarriers must be 1 or 2 (the state space grows as "
        "n^N/N!; use kinetic Monte Carlo for more carriers)");
  }
  ToFGraph graph(top, tof_, carriertype_, temperature_, field_);
  XTP_LOG(Log::error, log_) << graph.Summary() << std::flush;
  for (Index s : graph.DisconnectedSources()) {
    XTP_LOG(Log::info, log_)
        << "source segment " << s << " has no path to a sink, excluded"
        << std::flush;
  }

  GraphRandomWalk grw(graph);
  grw.setSettings(settings_);
  GraphRandomWalk::Result result = (carriers_ == 1)
                                       ? grw.SingleCarrier(solver_)
                                       : grw.TwoCarriers(observable_, solver_);
  bool transit =
      carriers_ == 2 && observable_ == GraphRandomWalk::Observable::Transit;
  if (carriers_ == 2) {
    XTP_LOG(Log::error, log_)
        << "Two carriers with site exclusion, observable: "
        << (transit ? "mean transit time per carrier (carriers removed at "
                      "the sink)"
                    : "time until the first carrier reaches the sink")
        << std::flush;
  }

  XTP_LOG(Log::error, log_)
      << "Solver: " << result.solver << ", " << result.unknowns << " unknowns, "
      << result.iterations << " iterations, backward error "
      << result.backward_error << " after " << result.refinement_steps
      << " refinement steps" << std::flush;
  if (!result.converged) {
    XTP_LOG(Log::error, log_)
        << "Warning: not converged (iteration limit or backward error above "
           "1e-10)"
        << std::flush;
  }
  XTP_LOG(Log::error, log_)
      << std::setprecision(6) << "Time of flight (uniform injection over "
      << result.sources.size() << " source "
      << (carriers_ == 1 ? "sites" : "pairs") << "): mean " << result.tau_mean
      << " s" << std::flush;
  if (!transit) {
    XTP_LOG(Log::error, log_) << "Standard deviation of the hitting time: "
                              << result.tau_std_injection << " s" << std::flush;
  }
  XTP_LOG(Log::error, log_)
      << "Harmonic mean over the injection states: " << result.tau_harmonic
      << " s" << std::flush;

  double Fn = field_Vm_.dot(graph.Normal());  // V/m along the normal
  double L = graph.Length() * tools::conv::bohr2nm * 1e-9;  // m
  if (Fn != 0.0) {
    double charge = 0.0;
    if (carriertype_ == QMStateType::Hole) {
      charge = 1.0;
    } else if (carriertype_ == QMStateType::Electron) {
      charge = -1.0;
    }
    if (charge * Fn <= 0.0) {
      XTP_LOG(Log::error, log_)
          << "Warning: the field does not drive the carriers towards the sink"
          << std::flush;
    }
    double factor = 1e4 / std::abs(Fn);  // m^2/Vs -> cm^2/Vs
    XTP_LOG(Log::error, log_)
        << "Mobility L/(tau F): " << L / result.tau_mean * factor
        << " cm^2/Vs (mean), " << L / result.tau_harmonic * factor
        << " cm^2/Vs (harmonic), with L = " << L * 1e9 << " nm and F = " << Fn
        << " V/m along the normal" << std::flush;
  }

  if (!resultfile_.empty()) {
    ResultFile out("grw: mean hitting times of the sink (graph random walk)");
    out.Add("carriers", result.carriers);
    out.Add("observable",
            transit ? std::string("transit") : std::string("first"));
    out.Add("unknowns", result.unknowns);
    out.Add("sources", Index(result.sources.size()));
    out.Add("tau_mean", result.tau_mean, "s");
    out.Add("tau_harmonic", result.tau_harmonic, "s");
    if (!transit) {
      out.Add("tau_std", result.tau_std_injection, "s");
    }
    out.Add("length", L, "m");
    out.Add("field_normal", Fn, "V/m");
    if (Fn != 0.0) {
      out.Add("mobility_mean", L / result.tau_mean * 1e4 / std::abs(Fn),
              "cm^2/Vs");
      out.Add("mobility_harmonic", L / result.tau_harmonic * 1e4 / std::abs(Fn),
              "cm^2/Vs");
    }
    out.Add("backward_error", result.backward_error, "-");
    out.Add("converged", Index(result.converged));
    out.Write(resultfile_);
  }
  if (!sourcefile_.empty()) {
    WriteSources(graph, result);
  }
  if (!occupationfile_.empty()) {
    WriteOccupation(graph, result);
  }
  std::cout << log_;
  return true;
}

void GRW::WriteSources(const ToFGraph& graph,
                       const GraphRandomWalk::Result& result) const {
  std::ofstream out(sourcefile_);
  out << "# per injection state: for each carrier segment id, x along the "
         "normal [nm], site energy [eV]; then the mean time [s] and, if "
         "available, the standard deviation of the hitting time [s]\n";
  out << std::setprecision(10);
  for (Index k = 0; k < Index(result.sources.size()); ++k) {
    for (Index s : result.sources[k]) {
      const ToFGraph::Site& site = graph.getSite(s);
      out << site.segment << " " << site.x * tools::conv::bohr2nm << " "
          << site.energy * tools::conv::hrt2ev << " ";
    }
    out << result.tau[k];
    if (result.tau_std.size() == result.tau.size()) {
      out << " " << result.tau_std[k];
    }
    out << "\n";
  }
}

void GRW::WriteOccupation(const ToFGraph& graph,
                          const GraphRandomWalk::Result& result) const {
  std::ofstream out(occupationfile_);
  out << "# expected time each segment is occupied before the carriers reach "
         "the sink, uniform injection\n";
  out << "# segment id, x along the normal [nm], site energy [eV], role, "
         "time [s], fraction of the total occupation time\n";
  out << std::setprecision(10);
  const char* roles[] = {"bulk", "source", "sink"};
  for (const ToFGraph::Site& site : graph.Sites()) {
    double t = result.occupation_time[site.segment];
    out << site.segment << " " << site.x * tools::conv::bohr2nm << " "
        << site.energy * tools::conv::hrt2ev << " "
        << roles[static_cast<int>(site.role)] << " " << t << " "
        << t / (double(result.carriers) * result.tau_mean) << "\n";
  }
}

}  // namespace xtp
}  // namespace votca
