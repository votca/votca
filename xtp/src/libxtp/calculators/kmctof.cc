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
#include <chrono>
#include <fstream>
#include <iomanip>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/eigen.h"
#include "votca/xtp/kmcrates.h"
#include "votca/xtp/kmcwaitingtime.h"
#include "votca/xtp/resultfile.h"
#include "votca/xtp/tofsimulation.h"
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "kmctof.h"

namespace votca {
namespace xtp {

void KMCToF::ParseOptions(const tools::Property& options) {
  carriertype_ = QMStateType(options.get(".carriertype").as<std::string>());
  if (!carriertype_.isKMCState()) {
    throw std::runtime_error("kmctof cannot be run for state: " +
                             carriertype_.ToLongString());
  }
  temperature_ = options.get(".temperature").as<double>() * tools::conv::kB *
                 tools::conv::ev2hrt;
  field_Vm_ = options.get(".field").as<Eigen::Vector3d>();
  field_ = field_Vm_ * (tools::conv::ev2hrt / (1e9 * tools::conv::nm2bohr));
  tof_ = ToFGraph::OptionsFromProperty(options.get(".tof"));
  seed_ = options.get(".seed").as<std::uint64_t>();
  runs_ = options.get(".runs").as<Index>();

  const tools::Property& inj = options.get(".injection");
  KMCEngine::Pulse pulse;
  pulse.carriers = inj.get("carriers").as<Index>();
  pulse.shape = KMCEngine::ShapeFromString(inj.get("shape").as<std::string>());
  pulse.width = inj.get("width").as<double>();
  Index npulses = inj.get("pulses").as<Index>();
  double period = inj.get("period").as<double>();
  if (npulses < 1 || pulse.carriers < 1) {
    throw std::runtime_error("kmctof: pulses and carriers must be positive");
  }
  if (npulses > 1 && !(period > 0.0)) {
    throw std::runtime_error("kmctof: a pulse train needs period > 0");
  }
  pulses_.clear();
  for (Index p = 0; p < npulses; ++p) {
    pulse.start = double(p) * period;
    pulses_.push_back(pulse);
  }

  waitingtime_ = options.get(".waitingtime").as<std::string>();
  limits_.max_time = options.get(".maxtime").as<double>();
  limits_.max_events = options.get(".maxevents").as<Index>();
  runfile_ = options.get(".runfile").as<std::string>();
  currentfile_ = options.get(".current.file").as<std::string>();
  resultfile_ = options.get(".resultfile").as<std::string>();
  tmin_ = options.get(".current.tmin").as<double>();
  tmax_ = options.get(".current.tmax").as<double>();
  bins_per_decade_ = options.get(".current.binsperdecade").as<Index>();
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n ...");
}

bool KMCToF::Evaluate(Topology& top) {
  XTP_LOG(Log::error, log_) << "\n-----------------------------------"
                               "\n  KINETIC MONTE CARLO TIME OF FLIGHT"
                               "\n-----------------------------------"
                            << std::flush;
  ToFGraph graph(top, tof_, carriertype_, temperature_, field_);
  XTP_LOG(Log::error, log_) << graph.Summary() << std::flush;

  ToFSimulation::Settings settings;
  settings.seed = seed_;
  settings.runs = runs_;
  settings.pulses = pulses_;
  settings.limits = limits_;
  settings.tmin = tmin_;
  settings.tmax = tmax_;
  settings.bins_per_decade = bins_per_decade_;
  // rates and waiting times: the static Marcus rates of the graph; a
  // time-dependent TimeDependentRates (e.g. oscillating field) plugs in here
  StaticRates rates(graph);
  std::unique_ptr<WaitingTimeStrategy> strategy =
      MakeWaitingTimeStrategy(waitingtime_);
  ToFSimulation simulation(graph, settings, rates, *strategy);
  XTP_LOG(Log::error, log_)
      << "Waiting times: " << simulation.StrategyName() << std::flush;
  XTP_LOG(Log::error, log_)
      << runs_ << " runs, " << pulses_.size() << " pulse(s) of "
      << pulses_[0].carriers << " carrier(s), " << OPENMP::getMaxThreads()
      << " thread(s)" << std::flush;
  auto start = std::chrono::steady_clock::now();
  ToFSimulation::Result result = simulation.Simulate();
  double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count();

  XTP_LOG(Log::error, log_)
      << std::setprecision(6) << result.hops << " hops, " << result.null_events
      << " null events (occupied site or rejected), " << seconds << " s"
      << std::flush;
  if (result.censored > 0) {
    XTP_LOG(Log::error, log_)
        << "Warning: " << result.censored
        << " runs stopped by maxtime/maxevents before all carriers arrived; "
           "they are left out, which biases the times low"
        << std::flush;
  }
  if (result.blocked_injections > 0) {
    XTP_LOG(Log::error, log_)
        << "Warning: " << result.blocked_injections
        << " injections dropped (all source sites occupied)" << std::flush;
  }
  auto report = [this](const std::string& name,
                       const ToFSimulation::Statistic& s) {
    XTP_LOG(Log::error, log_)
        << name << ": " << s.mean << " +- " << s.standard_error
        << " s (standard error, " << s.samples << " runs), standard deviation "
        << s.standard_deviation
        << " s, largest 1 % of the runs: " << 100.0 * s.top_percent_share
        << " % of the sum" << std::flush;
  };
  report("First arrival", result.first);
  report("Mean transit time per carrier", result.transit);
  if (result.transit.top_percent_share > 0.2) {
    XTP_LOG(Log::error, log_)
        << "Warning: heavy tail, the mean is dominated by rare long runs; "
           "compare with grw"
        << std::flush;
  }

  double Fn = field_Vm_.dot(graph.Normal());
  if (Fn != 0.0 && result.transit.samples > 0) {
    double L = graph.Length() * tools::conv::bohr2nm * 1e-9;  // m
    XTP_LOG(Log::error, log_)
        << "Mobility L/(tau F) from the mean transit time: "
        << L / result.transit.mean * 1e4 / std::abs(Fn) << " cm^2/Vs"
        << std::flush;
  }

  if (!resultfile_.empty()) {
    ResultFile out("kmctof: time-of-flight kinetic Monte Carlo");
    out.Add("runs", runs_);
    out.Add("injected", result.injected);
    out.Add("censored", result.censored);
    out.Add("blocked_injections", result.blocked_injections);
    out.Add("hops", result.hops);
    out.Add("null_events", result.null_events);
    out.Add("first_mean", result.first.mean, "s");
    out.Add("first_error", result.first.standard_error, "s");
    out.Add("first_std", result.first.standard_deviation, "s");
    out.Add("transit_mean", result.transit.mean, "s");
    out.Add("transit_error", result.transit.standard_error, "s");
    out.Add("transit_std", result.transit.standard_deviation, "s");
    out.Add("transit_top_percent_share", result.transit.top_percent_share, "-");
    double L = graph.Length() * tools::conv::bohr2nm * 1e-9;  // m
    out.Add("length", L, "m");
    out.Add("field_normal", Fn, "V/m");
    if (Fn != 0.0 && result.transit.samples > 0) {
      out.Add("mobility_transit", L / result.transit.mean * 1e4 / std::abs(Fn),
              "cm^2/Vs");
    }
    out.Write(resultfile_);
  }
  if (!runfile_.empty()) {
    std::ofstream out(runfile_);
    out << "# run, injected, absorbed, censored, first arrival [s], mean "
           "transit [s], hops, null events, end time [s]\n";
    out << std::setprecision(10);
    for (Index r = 0; r < Index(result.runs.size()); ++r) {
      const ToFSimulation::RunData& d = result.runs[r];
      out << r << " " << d.status.injected << " " << d.status.absorbed << " "
          << d.status.censored << " " << d.first_arrival << " "
          << d.mean_transit << " " << d.status.hops << " "
          << d.status.null_events << " " << d.status.end_time << "\n";
    }
  }
  if (!currentfile_.empty()) {
    std::ofstream out(currentfile_);
    out << "# transient current per injected carrier (Ramo-Shockley): "
           "I = (sum of the hop displacements along the transport normal "
           "in the bin) / (L x injected carriers x bin width), in 1/s; the "
           "circuit current is q x injected carriers x I, and the integral "
           "of I is the mean displacement per injected carrier in units of "
           "L\n";
    out << "# collected fraction: carriers absorbed up to t_high / injected "
           "carriers\n";
    out << "# t_low [s], t_high [s], t_mid (geometric) [s], current [1/s], "
           "collected fraction at t_high\n";
    out << "# displacement before t_low of the first bin: "
        << result.displacement_before
        << ", after the last bin: " << result.displacement_after << "\n";
    out << std::setprecision(10);
    for (Index b = 0; b < Index(result.current.size()); ++b) {
      double lo = result.bin_edges[b];
      double hi = result.bin_edges[b + 1];
      out << lo << " " << hi << " " << std::sqrt(lo * hi) << " "
          << result.current[b] << " " << result.collected[b] << "\n";
    }
  }
  std::cout << log_;
  return true;
}

}  // namespace xtp
}  // namespace votca
