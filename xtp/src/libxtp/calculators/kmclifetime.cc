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
#include <cmath>
#include <fstream>
#include <iomanip>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/eigen.h"
#include "votca/xtp/kmcrates.h"
#include "votca/xtp/kmcwaitingtime.h"
#include "votca/xtp/rate_engine.h"
#include "votca/xtp/resultfile.h"
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "kmclifetime.h"

namespace votca {
namespace xtp {

void KMCLifetime::ParseOptions(const tools::Property& options) {
  carriertype_ = QMStateType(options.get(".carriertype").as<std::string>());
  if (!carriertype_.isKMCState()) {
    throw std::runtime_error("kmclifetime cannot be run for state: " +
                             carriertype_.ToLongString());
  }
  temperature_ = options.get(".temperature").as<double>() * tools::conv::kB *
                 tools::conv::ev2hrt;
  injection_pattern_ = options.get(".injectionpattern").as<std::string>();
  ignore_types_ = options.get(".ignoresegments").as<std::string>();
  if (options.get(".injectionmethod").as<std::string>() != "random") {
    throw std::runtime_error("kmclifetime: injectionmethod must be random");
  }
  lifetimefile_ = options.get(".lifetimefile").as<std::string>();
  insertions_ = options.get(".numberofinsertions").as<Index>();
  settings_.seed = options.get(".seed").as<std::uint64_t>();
  settings_.carriers = options.get(".numberofcarriers").as<Index>();
  settings_.runs = options.get(".runs").as<Index>();
  settings_.decays = insertions_;
  settings_.max_wall_time = options.get(".maxrealtime").as<double>();
  exact_ = options.get(".exact").as<bool>();
  waitingtime_ = options.get(".waitingtime").as<std::string>();
  trajectoryfile_ = options.get(".trajectoryfile").as<std::string>();
  settings_.keep_records = !trajectoryfile_.empty();
  ratefile_ = options.get(".ratefile").as<std::string>();
  occfile_ = options.get(".occfile").as<std::string>();
  probfile_ = options.get(".decayprobfile").as<std::string>();
  resultfile_ = options.get(".resultfile").as<std::string>();

  const tools::Property& energy = options.get(".carrierenergy");
  carrierenergy_ = energy.get("run").as<bool>();
  energyfile_ = energy.get("outputfile").as<std::string>();
  alpha_ = energy.get("alpha").as<double>();
  if (carrierenergy_) {
    settings_.energy_hops = energy.get("outputsteps").as<Index>();
  }
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n ...");
}

std::vector<double> KMCLifetime::ReadDecayRates(const ToFGraph& graph) const {
  tools::Property xml;
  xml.LoadFromXML(lifetimefile_);
  std::vector<tools::Property*> sites = xml.Select("lifetimes.site");
  const Index n = graph.NumberOfSites();
  std::vector<double> rates(n, -1.0);
  for (tools::Property* prop : sites) {
    Index id = prop->getAttribute<Index>("id");
    if (id < 0 || id >= n) {
      throw std::runtime_error("kmclifetime: site " + std::to_string(id) +
                               " in " + lifetimefile_ +
                               " is not in the state file");
    }
    if (rates[id] >= 0.0) {
      throw std::runtime_error("kmclifetime: site " + std::to_string(id) +
                               " appears twice in " + lifetimefile_);
    }
    double lifetime = prop->as<double>();
    if (!(lifetime > 0.0) || !std::isfinite(lifetime)) {
      throw std::runtime_error("kmclifetime: the lifetime of site " +
                               std::to_string(id) + " must be positive");
    }
    rates[id] = 1.0 / lifetime;
  }
  for (Index i = 0; i < n; ++i) {
    if (rates[i] < 0.0) {
      if (!graph.getSite(i).ignored) {
        throw std::runtime_error("kmclifetime: no lifetime for site " +
                                 std::to_string(i) + " in " + lifetimefile_);
      }
      rates[i] = 0.0;
    }
  }
  return rates;
}

bool KMCLifetime::Evaluate(Topology& top) {
  XTP_LOG(Log::error, log_) << "\n-----------------------------------"
                               "\n      KMC WITH FINITE LIFETIMES"
                               "\n-----------------------------------"
                            << std::flush;
  ToFGraph graph(top, ToFGraph::BulkOptions(injection_pattern_, ignore_types_),
                 carriertype_, temperature_, Eigen::Vector3d::Zero());
  XTP_LOG(Log::error, log_) << graph.Summary() << std::flush;
  std::vector<double> decay = ReadDecayRates(graph);
  if (!ratefile_.empty()) {
    WriteRates(top);
  }
  double nm = tools::conv::bohr2nm;

  std::unique_ptr<LifetimeSimulation::ExactResult> exact;
  if (exact_) {
    exact = std::make_unique<LifetimeSimulation::ExactResult>(
        LifetimeSimulation::Exact(graph, decay));
    XTP_LOG(Log::error, log_) << "Exact for one " << carriertype_.ToLongString()
                              << " (" << exact->solver << ", backward error "
                              << exact->backward_error << "):" << std::flush;
    XTP_LOG(Log::error, log_)
        << std::scientific << "  lifetime " << exact->lifetime
        << " s, diffusion length sqrt(<R^2>) "
        << std::sqrt(exact->squared_displacement) * nm
        << " nm, per direction sqrt(<R_a^2>) "
        << (exact->squared_components.cwiseSqrt() * nm).transpose() << " nm"
        << std::flush;
    if (!exact->converged) {
      XTP_LOG(Log::error, log_)
          << "Warning: the iterative solver did not converge" << std::flush;
    }
    if (settings_.carriers > 1) {
      XTP_LOG(Log::error, log_)
          << "  (with " << settings_.carriers
          << " carriers the KMC includes site exclusion, the exact values "
             "do not)"
          << std::flush;
    }
  }

  std::unique_ptr<LifetimeSimulation::Result> result;
  if (insertions_ > 0) {
    StaticRates rates(graph);
    std::unique_ptr<WaitingTimeStrategy> strategy =
        MakeWaitingTimeStrategy(waitingtime_);
    LifetimeSimulation simulation(graph, decay, settings_, rates, *strategy);
    XTP_LOG(Log::error, log_)
        << settings_.carriers << " " << carriertype_.ToLongString()
        << "(s) at a time, " << insertions_ << " decays in " << settings_.runs
        << " run(s), " << OPENMP::getMaxThreads() << " thread(s), waiting "
        << "times: " << strategy->Name() << std::flush;
    auto start = std::chrono::steady_clock::now();
    result =
        std::make_unique<LifetimeSimulation::Result>(simulation.Simulate());
    double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    const LifetimeSimulation::Result& r = *result;
    XTP_LOG(Log::error, log_)
        << std::setprecision(6) << r.decays << " decays, " << r.hop_events
        << " hops, " << r.null_events << " attempts into occupied sites, "
        << r.simulated_time << " s simulated, " << seconds << " s"
        << std::flush;
    if (r.wall_time_stops > 0) {
      XTP_LOG(Log::error, log_)
          << "Warning: " << r.wall_time_stops
          << " run(s) stopped by maxrealtime" << std::flush;
    }
    double LD = std::sqrt(r.squared_displacement);
    XTP_LOG(Log::error, log_)
        << std::scientific << "Average lifetime: " << r.lifetime << " +- "
        << r.lifetime_error << " s (standard error)\n"
        << "  hops per " << carriertype_.ToLongString() << ": " << r.hops
        << "\n  mean distance <|R|>: " << r.mean_distance * nm << " +- "
        << r.mean_distance_error * nm << " nm"
        << "\n  diffusion length sqrt(<R^2>): " << LD * nm << " +- "
        << r.squared_displacement_error / (2.0 * LD) * nm << " nm"
        << "\n  per direction sqrt(<R_a^2>): "
        << (r.squared_components.cwiseSqrt() * nm).transpose() << " nm"
        << "\n  D = <R^2>/(6 <lifetime>): "
        << r.squared_displacement / (6.0 * r.lifetime) * nm * nm << " nm^2/s"
        << std::flush;
    if (exact && settings_.carriers == 1) {
      XTP_LOG(Log::error, log_)
          << "  KMC - exact: lifetime "
          << (r.lifetime - exact->lifetime) / r.lifetime_error
          << " standard errors, <R^2> "
          << (r.squared_displacement - exact->squared_displacement) /
                 r.squared_displacement_error
          << " standard errors" << std::flush;
    }

    if (!trajectoryfile_.empty()) {
      std::ofstream out(trajectoryfile_);
      out << "#Simtime [s]\t Insertion\t Carrier ID\t Lifetime[s]\tSteps\t "
             "Decay segment\t x_travelled[nm]\t y_travelled[nm]\t "
             "z_travelled[nm]\t run\n";
      out << std::setprecision(10);
      for (Index i = 0; i < Index(r.records.size()); ++i) {
        const LifetimeSimulation::Record& d = r.records[i];
        out << d.time << "\t" << i << "\t" << d.carrier << "\t" << d.lifetime
            << "\t" << d.hops << "\t" << d.site << "\t"
            << d.displacement.x() * nm << "\t" << d.displacement.y() * nm
            << "\t" << d.displacement.z() * nm << "\t" << d.run << "\n";
      }
      XTP_LOG(Log::error, log_)
          << "Decays written to " << trajectoryfile_ << std::flush;
    }
    if (carrierenergy_) {
      std::ofstream out(energyfile_);
      out << "#Simtime [s]\tHops\tEnergy[eV]\tEnergy_a=" << alpha_
          << "[eV] (run 0, mean over the carriers present)\n";
      double avg = r.energy.empty() ? 0.0 : r.energy.front().energy;
      for (const LifetimeSimulation::EnergySample& e : r.energy) {
        avg = alpha_ * e.energy + (1.0 - alpha_) * avg;
        out << e.time << "\t" << e.hops << "\t"
            << e.energy * tools::conv::hrt2ev << "\t"
            << avg * tools::conv::hrt2ev << "\n";
      }
    }
    if (!occfile_.empty()) {
      std::ofstream out(occfile_);
      out << "#SiteID, Occupation prob at "
          << temperature_ * tools::conv::hrt2ev / tools::conv::kB
          << "K for carrier:" << carriertype_.ToString() << "\n";
      for (Index i = 0; i < graph.NumberOfSites(); ++i) {
        out << i << "\t" << r.occupation[i] << "\n";
      }
    }
  }
  if (!probfile_.empty()) {
    WriteDecayProbabilities(graph, decay, result.get(), exact.get());
  }
  if (!resultfile_.empty()) {
    ResultFile out("kmclifetime: kinetic Monte Carlo with finite lifetimes");
    const double to_m = nm * 1e-9;
    const char* axes = "xyz";
    out.Add("carriers", settings_.carriers);
    if (result) {
      const LifetimeSimulation::Result& r = *result;
      out.Add("decays", r.decays);
      out.Add("hop_events", r.hop_events);
      out.Add("null_events", r.null_events);
      out.Add("lifetime", r.lifetime, "s");
      out.Add("lifetime_error", r.lifetime_error, "s");
      out.Add("hops_per_carrier", r.hops, "-");
      out.Add("mean_distance", r.mean_distance * to_m, "m");
      out.Add("mean_distance_error", r.mean_distance_error * to_m, "m");
      out.Add("msd", r.squared_displacement * to_m * to_m, "m^2");
      out.Add("msd_error", r.squared_displacement_error * to_m * to_m, "m^2");
      for (Index a = 0; a < 3; ++a) {
        out.Add(std::string("msd_") + axes[a],
                r.squared_components[a] * to_m * to_m, "m^2");
      }
      out.Add("diffusion_length", std::sqrt(r.squared_displacement) * to_m,
              "m");
      out.Add("diffusion_coefficient",
              r.squared_displacement / (6.0 * r.lifetime) * to_m * to_m,
              "m^2/s");
    }
    out.Add("exact_available", Index(exact != nullptr));
    if (exact) {
      out.Add("exact_lifetime", exact->lifetime, "s");
      out.Add("exact_msd", exact->squared_displacement * to_m * to_m, "m^2");
      for (Index a = 0; a < 3; ++a) {
        out.Add(std::string("exact_msd_") + axes[a],
                exact->squared_components[a] * to_m * to_m, "m^2");
      }
      out.Add("exact_diffusion_length",
              std::sqrt(exact->squared_displacement) * to_m, "m");
      out.Add("exact_backward_error", exact->backward_error, "-");
    }
    out.Write(resultfile_);
  }
  std::cout << log_;
  return true;
}

void KMCLifetime::WriteRates(const Topology& top) {
  std::ofstream out(ratefile_);
  out << "#PairID,SiteID1,SiteID2, ,rate12[1/s],rate21[1/s] at "
      << temperature_ * tools::conv::hrt2ev / tools::conv::kB
      << "K for carrier:" << carriertype_.ToString() << "\n";
  Rate_Engine engine(temperature_, Eigen::Vector3d::Zero());
  for (const QMPair* pair : top.NBList()) {
    Rate_Engine::PairRates rates = engine.Rate(*pair, carriertype_);
    out << pair->getId() << " " << pair->Seg1()->getId() << " "
        << pair->Seg2()->getId() << " " << rates.rate12 << " " << rates.rate21
        << "\n";
  }
  XTP_LOG(Log::error, log_) << "Rates written to " << ratefile_ << std::flush;
}

void KMCLifetime::WriteDecayProbabilities(
    const ToFGraph& graph, const std::vector<double>& decay,
    const LifetimeSimulation::Result* kmc,
    const LifetimeSimulation::ExactResult* exact) {
  const Index n = graph.NumberOfSites();
  Eigen::VectorXd in = Eigen::VectorXd::Zero(n);
  for (Index i = 0; i < n; ++i) {
    for (const ToFGraph::Edge& e : graph.OutEdges(i)) {
      in[e.to] += e.rate;
    }
  }
  std::ofstream out(probfile_);
  out << "#SiteID, k/(k+outgoing rates), k/(k+incoming rates), fraction of "
         "the KMC decays, exact decay probability (one carrier)\n";
  for (Index i = 0; i < n; ++i) {
    double k = decay[i];
    double pout = (k > 0.0) ? k / (k + graph.EscapeRate(i)) : 0.0;
    double pin = (k > 0.0) ? k / (k + in[i]) : 0.0;
    out << i << " " << pout << " " << pin << " "
        << (kmc != nullptr ? kmc->decay_fraction[i] : 0.0) << " "
        << (exact != nullptr ? exact->decay_probability[i] : 0.0) << "\n";
  }
  XTP_LOG(Log::error, log_)
      << "Decay probabilities written to " << probfile_ << std::flush;
}

}  // namespace xtp
}  // namespace votca
