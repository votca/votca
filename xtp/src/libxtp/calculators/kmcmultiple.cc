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
#include "kmcmultiple.h"

namespace votca {
namespace xtp {

void KMCMultiple::ParseOptions(const tools::Property& options) {
  carriertype_ = QMStateType(options.get(".carriertype").as<std::string>());
  if (!carriertype_.isKMCState()) {
    throw std::runtime_error("kmcmultiple cannot be run for state: " +
                             carriertype_.ToLongString());
  }
  temperature_ = options.get(".temperature").as<double>() * tools::conv::kB *
                 tools::conv::ev2hrt;
  field_Vm_ = options.get(".field").as<Eigen::Vector3d>();
  field_ = field_Vm_ * (tools::conv::ev2hrt / (1e9 * tools::conv::nm2bohr));
  injection_pattern_ = options.get(".injectionpattern").as<std::string>();
  ignore_types_ = options.get(".ignoresegments").as<std::string>();
  std::string method = options.get(".injectionmethod").as<std::string>();
  if (method != "random") {
    throw std::runtime_error("kmcmultiple: injectionmethod must be random");
  }

  settings_.seed = options.get(".seed").as<std::uint64_t>();
  settings_.carriers = options.get(".numberofcarriers").as<Index>();
  settings_.runs = options.get(".runs").as<Index>();
  settings_.blocks = options.get(".blocks").as<Index>();
  settings_.equilibration = options.get(".equilibration").as<double>();
  settings_.limits.max_wall_time = options.get(".maxrealtime").as<double>();

  runtime_ = options.get(".runtime").as<double>();
  outputtime_ = options.get(".outputtime").as<double>();
  steps_ = runtime_ > 100.0;
  if (steps_) {
    if (outputtime_ != std::floor(outputtime_)) {
      throw std::runtime_error(
          "kmcmultiple: runtime is a number of events (> 100), so outputtime "
          "must be an integer number of hops too");
    }
    settings_.limits.max_events = Index(runtime_);
    settings_.sample_hops = Index(outputtime_);
  } else {
    if (!(runtime_ > 0.0)) {
      throw std::runtime_error("kmcmultiple: runtime must be positive");
    }
    settings_.limits.max_time = runtime_;
    settings_.sample_interval = outputtime_;
  }

  waitingtime_ = options.get(".waitingtime").as<std::string>();
  trajectoryfile_ = options.get(".trajectoryfile").as<std::string>();
  timefile_ = options.get(".timefile").as<std::string>();
  ratefile_ = options.get(".ratefile").as<std::string>();
  occfile_ = options.get(".occfile").as<std::string>();
  resultfile_ = options.get(".resultfile").as<std::string>();
  exact_ = options.get(".exact").as<bool>();
  if (outputtime_ == 0.0) {
    trajectoryfile_.clear();
    timefile_.clear();
  }
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n ...");
}

double KMCMultiple::Charge() const {
  if (carriertype_ == QMStateType::Electron) {
    return -1.0;
  } else if (carriertype_ == QMStateType::Hole) {
    return 1.0;
  }
  return 0.0;
}

bool KMCMultiple::Evaluate(Topology& top) {
  XTP_LOG(Log::error, log_) << "\n-----------------------------------"
                               "\n      KMC FOR MULTIPLE CHARGES"
                               "\n-----------------------------------"
                            << std::flush;
  ToFGraph::Options opt =
      ToFGraph::BulkOptions(injection_pattern_, ignore_types_);
  ToFGraph graph(top, opt, carriertype_, temperature_, field_);
  XTP_LOG(Log::error, log_) << graph.Summary() << std::flush;
  if (!ratefile_.empty()) {
    WriteRates(top);
  }

  double nm = tools::conv::bohr2nm;
  double volume = top.BoxVolume() * nm * nm * nm;
  XTP_LOG(Log::error, log_)
      << settings_.carriers << " " << carriertype_.ToLongString()
      << "(s), density " << double(settings_.carriers) / volume << " nm^-3, "
      << settings_.runs << " run(s), " << OPENMP::getMaxThreads()
      << " thread(s)" << std::flush;
  if (steps_) {
    XTP_LOG(Log::error, log_)
        << "Each run stops after " << settings_.limits.max_events
        << " KMC events (hops and attempts into occupied sites)" << std::flush;
  } else {
    XTP_LOG(Log::error, log_) << "Each run stops at " << runtime_
                              << " s simulated time" << std::flush;
  }

  StaticRates rates(graph);
  std::unique_ptr<WaitingTimeStrategy> strategy =
      MakeWaitingTimeStrategy(waitingtime_);
  BulkSimulation simulation(graph, settings_, rates, *strategy);
  XTP_LOG(Log::error, log_)
      << "Waiting times: " << strategy->Name() << std::flush;
  auto start = std::chrono::steady_clock::now();
  BulkSimulation::Result result = simulation.Simulate();
  double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count();

  double attempts = double(result.hops + result.null_events);
  XTP_LOG(Log::error, log_)
      << std::setprecision(6) << result.hops << " hops, " << result.null_events
      << " attempts into occupied sites ("
      << 100.0 * double(result.null_events) / std::max(attempts, 1.0) << " %), "
      << seconds << " s" << std::flush;
  if (result.wall_time_stops > 0) {
    XTP_LOG(Log::error, log_) << "Warning: " << result.wall_time_stops
                              << " run(s) stopped by maxrealtime" << std::flush;
  }
  const BulkSimulation::RunData& run0 = result.runs.front();
  XTP_LOG(Log::error, log_)
      << "Observation window of run 0: " << run0.start << " s to " << run0.end
      << " s, " << run0.blocks << " blocks of " << run0.block_length << " s"
      << std::flush;

  std::string error_kind = settings_.runs > 1
                               ? "standard error over runs"
                               : "standard error from the time blocks";
  XTP_LOG(Log::error, log_)
      << std::scientific
      << "Drift velocity (nm/s): " << (result.velocity * nm).transpose()
      << "\n    +- " << (result.velocity_error * nm).transpose() << " ("
      << error_kind << ")" << std::flush;

  // mobility along the field, mu = q v.F / |F|^2 (positive if the carriers
  // drift along the force); m^2/Vs
  const double F = field_Vm_.norm();
  const double to_m = nm * 1e-9;
  const bool driven = F > 0.0 && Charge() != 0.0;
  const Eigen::Vector3d dir =
      F > 0.0 ? Eigen::Vector3d(field_Vm_ / F) : Eigen::Vector3d::Zero();
  auto mobility = [&](const Eigen::Vector3d& v) {
    return Charge() * v.dot(dir) * to_m / F;
  };
  double mu = 0.0;
  double mu_err = 0.0;
  if (driven) {
    mu = mobility(result.velocity);
    mu_err =
        std::sqrt((result.velocity_error.cwiseProduct(dir)).squaredNorm()) *
        to_m / F;
    XTP_LOG(Log::error, log_)
        << "Mobility along the field: " << mu * 1e4 << " +- " << mu_err * 1e4
        << " cm^2/Vs (" << mu * 1e18 << " nm^2/Vs)" << std::flush;
  }

  Eigen::Matrix3d D = result.diffusion * nm * nm;
  XTP_LOG(Log::error, log_) << "Diffusion tensor (nm^2/s):\n"
                            << D << std::flush;
  if (settings_.runs > 1) {
    XTP_LOG(Log::error, log_) << "  +- (standard error over runs)\n"
                              << result.diffusion_error * nm * nm << std::flush;
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(D);
  for (Index i = 0; i < 3; ++i) {
    XTP_LOG(Log::error, log_)
        << "  eigenvalue " << es.eigenvalues()(i) << " nm^2/s, direction "
        << es.eigenvectors().col(i).transpose() << std::flush;
  }
  const double kT = temperature_ * tools::conv::hrt2ev;  // eV
  auto einstein = [kT](const Eigen::Matrix3d& Dnm) {
    return Dnm.trace() / 3.0 * 1e-18 / kT;  // m^2/Vs
  };
  if (F == 0.0) {
    XTP_LOG(Log::error, log_)
        << "Einstein mobility D/kT from the mean eigenvalue (isotropic "
           "medium, nondegenerate carriers): "
        << einstein(D) * 1e4 << " cm^2/Vs (" << einstein(D) * 1e18
        << " nm^2/Vs)" << std::flush;
  }

  BulkSimulation::ExactResult exact;
  if (exact_) {
    exact = BulkSimulation::Exact(graph);
    if (!exact.available) {
      XTP_LOG(Log::error, log_)
          << "Exact single-carrier values not available: " << exact.message
          << std::flush;
    } else {
      Eigen::Matrix3d Dx = exact.diffusion * nm * nm;
      XTP_LOG(Log::error, log_)
          << "Exact for one carrier (master equation, " << exact.solver << ", "
          << exact.sites << " sites, backward error " << exact.backward_error
          << "):\n  drift velocity (nm/s): "
          << (exact.velocity * nm).transpose() << "\n  diffusion tensor "
          << "(nm^2/s):\n"
          << Dx << std::flush;
      if (driven) {
        XTP_LOG(Log::error, log_)
            << "  mobility along the field: " << mobility(exact.velocity) * 1e4
            << " cm^2/Vs" << std::flush;
      } else if (F == 0.0) {
        XTP_LOG(Log::error, log_)
            << "  Einstein mobility: " << einstein(Dx) * 1e4 << " cm^2/Vs"
            << std::flush;
      }
      if (!exact.converged) {
        XTP_LOG(Log::error, log_)
            << "  Warning: the linear solver did not converge" << std::flush;
      }
      if (settings_.carriers > 1) {
        XTP_LOG(Log::error, log_)
            << "  (with " << settings_.carriers
            << " carriers the KMC includes site exclusion, the exact values "
               "are for an isolated carrier)"
            << std::flush;
      } else if (driven && mu_err > 0.0) {
        XTP_LOG(Log::error, log_) << "  KMC - exact mobility: "
                                  << (mu - mobility(exact.velocity)) / mu_err
                                  << " standard errors" << std::flush;
      }
    }
  }

  if (!resultfile_.empty()) {
    ResultFile out("kmcmultiple: bulk kinetic Monte Carlo");
    out.Add("carriers", settings_.carriers);
    out.Add("runs", settings_.runs);
    out.Add("hops", result.hops);
    out.Add("null_events", result.null_events);
    out.Add("window_start", run0.start, "s");
    out.Add("window_end", run0.end, "s");
    const char* axes = "xyz";
    for (Index a = 0; a < 3; ++a) {
      out.Add(std::string("velocity_") + axes[a], result.velocity[a] * to_m,
              "m/s");
      out.Add(std::string("velocity_error_") + axes[a],
              result.velocity_error[a] * to_m, "m/s");
    }
    if (driven) {
      out.Add("mobility", mu * 1e4, "cm^2/Vs");
      out.Add("mobility_error", mu_err * 1e4, "cm^2/Vs");
    }
    auto tensor = [&out, axes](const std::string& name,
                               const Eigen::Matrix3d& T) {
      for (Index a = 0; a < 3; ++a) {
        for (Index b = a; b < 3; ++b) {
          out.Add(name + "_" + axes[a] + axes[b], T(a, b) * 1e-18, "m^2/s");
        }
      }
    };
    tensor("D", D);
    if (settings_.runs > 1) {
      tensor("D_error", result.diffusion_error * nm * nm);
    }
    if (F == 0.0) {
      out.Add("einstein_mobility", einstein(D) * 1e4, "cm^2/Vs");
    }
    out.Add("exact_available", Index(exact.available));
    if (exact.available) {
      for (Index a = 0; a < 3; ++a) {
        out.Add(std::string("exact_velocity_") + axes[a],
                exact.velocity[a] * to_m, "m/s");
      }
      if (driven) {
        out.Add("exact_mobility", mobility(exact.velocity) * 1e4, "cm^2/Vs");
      }
      Eigen::Matrix3d Dx = exact.diffusion * nm * nm;
      tensor("exact_D", Dx);
      if (F == 0.0) {
        out.Add("exact_einstein_mobility", einstein(Dx) * 1e4, "cm^2/Vs");
      }
    }
    out.Write(resultfile_);
  }

  if (!occfile_.empty()) {
    WriteOccupation(graph, result);
  }
  if (!trajectoryfile_.empty()) {
    WriteTrajectory(result);
  }
  if (!timefile_.empty()) {
    WriteTimeFile(result);
  }
  std::cout << log_;
  return true;
}

void KMCMultiple::WriteRates(const Topology& top) {
  std::ofstream out(ratefile_);
  out << "#PairID,SiteID1,SiteID2, ,rate12[1/s],rate21[1/s] at "
      << temperature_ * tools::conv::hrt2ev / tools::conv::kB
      << "K for carrier:" << carriertype_.ToString() << "\n";
  Rate_Engine engine(temperature_, field_);
  for (const QMPair* pair : top.NBList()) {
    Rate_Engine::PairRates rates = engine.Rate(*pair, carriertype_);
    out << pair->getId() << " " << pair->Seg1()->getId() << " "
        << pair->Seg2()->getId() << " " << rates.rate12 << " " << rates.rate21
        << "\n";
  }
  XTP_LOG(Log::error, log_) << "Rates written to " << ratefile_ << std::flush;
}

void KMCMultiple::WriteOccupation(const ToFGraph& graph,
                                  const BulkSimulation::Result& result) {
  std::ofstream out(occfile_);
  out << "#SiteID, Occupation prob at "
      << temperature_ * tools::conv::hrt2ev / tools::conv::kB
      << "K for carrier:" << carriertype_.ToString() << "\n";
  for (Index i = 0; i < graph.NumberOfSites(); ++i) {
    out << i << "\t" << result.occupation[i] << "\n";
  }
  XTP_LOG(Log::error, log_)
      << "Occupation probabilities written to " << occfile_ << std::flush;
}

void KMCMultiple::WriteTrajectory(const BulkSimulation::Result& result) {
  std::ofstream out(trajectoryfile_);
  out << "time[s]\tsteps";
  for (Index c = 0; c < settings_.carriers; ++c) {
    out << "\tcarrier" << c + 1 << " x\tcarrier" << c + 1 << " y\tcarrier"
        << c + 1 << " z";
  }
  out << "\n";
  double nm = tools::conv::bohr2nm;
  for (const BulkSimulation::Sample& s : result.trajectory) {
    out << s.time << "\t" << s.hops;
    for (const Eigen::Vector3d& r : s.positions) {
      out << "\t" << r.x() * nm << "\t" << r.y() * nm << "\t" << r.z() * nm;
    }
    out << "\n";
  }
  XTP_LOG(Log::error, log_)
      << "Trajectory of run 0 written to " << trajectoryfile_ << std::flush;
}

void KMCMultiple::WriteTimeFile(const BulkSimulation::Result& result) {
  std::ofstream out(timefile_);
  out << "time[s]\tsteps\tenergy_per_carrier[eV]\tmobility[nm**2/Vs]\t"
         "distance_fielddirection[nm]\tdistance_absolute[nm]\n";
  if (result.trajectory.empty()) {
    return;
  }
  double nm = tools::conv::bohr2nm;
  double F = field_Vm_.norm();
  const std::vector<Eigen::Vector3d>& r0 = result.trajectory.front().positions;
  double t0 = result.trajectory.front().time;
  for (const BulkSimulation::Sample& s : result.trajectory) {
    Eigen::Vector3d dr = Eigen::Vector3d::Zero();
    for (Index c = 0; c < Index(s.positions.size()); ++c) {
      dr += s.positions[c] - r0[c];
    }
    dr *= nm / double(s.positions.size());  // nm, mean per carrier
    double along = (F > 0.0) ? dr.dot(field_Vm_) / F : 0.0;
    double mobility = 0.0;
    if (F > 0.0 && s.time > t0) {
      // nm/s over V/nm
      mobility = Charge() * along / (s.time - t0) / (F * 1e-9);
    }
    out << s.time << "\t" << s.hops << "\t" << s.energy * tools::conv::hrt2ev
        << "\t" << mobility << "\t" << along << "\t" << dr.norm() << "\n";
  }
}

}  // namespace xtp
}  // namespace votca
