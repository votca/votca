/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
#include <chrono>
#include <sstream>

// Third party includes
#include <boost/filesystem.hpp>
#include <numeric>
#include <stdexcept>

// Local VOTCA includes
#include "votca/tools/property.h"
#include "votca/tools/tokenizer.h"
#include "votca/xtp/jobtopology.h"
#include "votca/xtp/qmregion.h"

// Local private VOTCA includes
#include "qmmm.h"
#include "votca/xtp/qmstate.h"

namespace votca {
namespace xtp {

void QMMM::ParseSpecificOptions(const tools::Property& options) {

  print_regions_pdb_ = options.get(".print_regions_pdb").as<bool>();
  max_iterations_ = options.get(".max_iterations").as<Index>();
  regions_def_.second = options.get(".regions");
  regions_def_.first = mapfile_;
  use_gs_for_ex_ = options.get(".use_gs_for_ex").as<bool>();

  states_ = options.get(".io_states").as<std::vector<QMState>>();
  whichSegments_ = options.get(".segments").as<std::string>();

  if (whichSegments_ == "all") {
    all_segments_ = true;
  } else {
    all_segments_ = false;
    std::stringstream ss(whichSegments_);
    Index segID;
    while (ss >> segID) {
      segments_.push_back(segID);
    }
  }

  bool groundstate_found = std::any_of(
      states_.begin(), states_.end(),
      [](const QMState& state) { return state.Type() == QMStateType::Gstate; });
  if (!groundstate_found) {
    states_.push_back(QMState("n"));
  }
}

Job::JobResult QMMM::EvalJob(const Topology& top, Job& job, QMThread& Thread) {
  std::chrono::time_point<std::chrono::system_clock> start =
      std::chrono::system_clock::now();

  std::string qmmm_work_dir = "QMMM";
  if (!this->hasQMRegion()) {
    qmmm_work_dir = "MMMM";
  }
  std::string frame_dir =
      "frame_" + boost::lexical_cast<std::string>(top.getStep());
  std::string job_dir =
      "job_" + std::to_string(job.getId()) + "_" + job.getTag();
  boost::filesystem::path arg_path;
  std::string workdir =
      (arg_path / qmmm_work_dir / frame_dir / job_dir).generic_string();
  boost::filesystem::create_directories(workdir);
  Job::JobResult jres = Job::JobResult();
  Logger& pLog = Thread.getLogger();

  JobTopology jobtop = JobTopology(job, pLog, workdir);
  if (job.getInput().exists("restart")) {
    std::string checkptfile = job.getInput().get("restart").as<std::string>();
    XTP_LOG(Log::error, pLog)
        << TimeStamp() << " Restart job from " << checkptfile << std::flush;
    jobtop.ReadFromHdf5(checkptfile);
  } else {
    jobtop.BuildRegions(top, regions_def_);
  }

  if (print_regions_pdb_) {
    std::string pdb_filename = "regions.pdb";
    XTP_LOG(Log::error, pLog) << TimeStamp() << " Writing jobtopology to "
                              << (workdir + "/" + pdb_filename) << std::flush;
    jobtop.WriteToPdb(workdir + "/" + pdb_filename);
  }

  Index no_static_regions = std::accumulate(
      jobtop.begin(), jobtop.end(), 0, [](Index count, const auto& region) {
        return count += Index(region->Converged());
      });

  bool no_top_scf = false;
  if (jobtop.size() - no_static_regions < 2) {
    XTP_LOG(Log::error, pLog)
        << TimeStamp() << " Only " << jobtop.size() - no_static_regions
        << " scf region is used. The remaining regions are static. So no "
           "inter regions scf is required. "
        << std::flush;
    no_top_scf = true;
    max_iterations_ = 1;
  }
  Index iteration = 0;
  for (; iteration < max_iterations_; iteration++) {

    XTP_LOG(Log::error, pLog)
        << TimeStamp() << " --Inter Region SCF Iteration " << iteration + 1
        << " of " << max_iterations_ << std::flush;

    for (std::unique_ptr<Region>& region : jobtop) {
      XTP_LOG(Log::error, pLog)
          << TimeStamp() << " Evaluating " << region->identify() << " "
          << region->getId() << std::flush;
      region->Reset();
      region->Evaluate(jobtop.Regions());
      if (!region->Successful()) {
        jres.setStatus(Job::JobStatus::FAILED);
        jres.setError(region->ErrorMsg());
        return jres;
      }
    }

    std::string checkpointfilename =
        "checkpoint_iter_" + std::to_string(iteration + 1) + ".hdf5";
    XTP_LOG(Log::error, pLog) << TimeStamp() << " Writing checkpoint to "
                              << checkpointfilename << std::flush;
    jobtop.WriteToHdf5(workdir + "/" + checkpointfilename);

    if (!no_top_scf) {
      std::vector<bool> converged_regions;
      for (std::unique_ptr<Region>& region : jobtop) {
        converged_regions.push_back(region->Converged());
      }

      double etot = std::accumulate(
          jobtop.begin(), jobtop.end(), 0.0,
          [](double e, const auto& region) { return region->Etotal() + e; });

      XTP_LOG(Log::error, pLog) << TimeStamp() << " --Total Energy all regions "
                                << etot << std::flush;

      bool all_regions_converged =
          std::all_of(converged_regions.begin(), converged_regions.end(),
                      [](bool i) { return i; });

      if (all_regions_converged) {
        XTP_LOG(Log::error, pLog)
            << TimeStamp() << " Job converged after " << iteration + 1
            << " iterations." << std::flush;
        jres.setStatus(Job::JobStatus::COMPLETE);
        break;
      }
      if (iteration == max_iterations_ - 1) {
        XTP_LOG(Log::error, pLog)
            << TimeStamp() << " Job did not converge after " << iteration + 1
            << " iterations.\n Writing results to jobfile." << std::flush;
        jres.setStatus(Job::JobStatus::FAILED);
        jres.setError("Inter Region SCF did not converge in " +
                      std::to_string(max_iterations_) + " iterations.");
      }
    } else {
      jres.setStatus(Job::JobStatus::COMPLETE);
    }
  }

  tools::Property results;
  tools::Property& jobresult = results.add("output", "");
  tools::Property& regionsresults = jobresult.add("regions", "");
  double etot = 0.0;
  double charge = 0.0;
  for (const std::unique_ptr<Region>& reg : jobtop) {
    reg->AddResults(regionsresults);
    etot += reg->Etotal();
    charge += reg->charge();
  }
  std::chrono::time_point<std::chrono::system_clock> end =
      std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;
  jobresult.add("E_tot", std::to_string(etot * tools::conv::hrt2ev));
  jobresult.add("Compute_Time", std::to_string(Index(elapsed_time.count())));
  jobresult.add("Total_Charge", std::to_string(charge));
  if (!no_top_scf) {
    jobresult.add("Iterations", std::to_string(iteration + 1));
  }
  jres.setOutput(results);
  return jres;
}

bool QMMM::hasQMRegion() const {
  Logger log;
  QMRegion QMdummy(0, log, "");
  return std::any_of(regions_def_.second.begin(), regions_def_.second.end(),
                     [&](const tools::Property& reg) {
                       return reg.name() == QMdummy.identify();
                     });
}

std::string QMMM::getFirstRegionName() const {
  for (const auto& reg : regions_def_.second) {
    if (reg.get("id").as<Index>() == 0) {
      return reg.name();
    }
  }
  throw std::runtime_error("region ids do not start at 0");
  return "";
}

void QMMM::WriteJobFile(const Topology& top) {

  std::cout << std::endl
            << "... ... Writing job file " << jobfile_ << std::flush;

  std::ofstream ofs;
  ofs.open(jobfile_, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("\nERROR: bad file handle: " + jobfile_);
  }

  ofs << "<jobs>" << std::endl;
  Index jobid = 0;
  if (all_segments_) {
    for (const Segment& seg : top.Segments()) {
      for (const QMState& state : states_) {
        Job job = createJob(seg, state, jobid);
        job.ToStream(ofs);
        jobid++;
      }
    }
  } else {
    for (Index segID : segments_) {
      const Segment& seg = top.Segments()[segID];
      for (const QMState& state : states_) {
        Job job = createJob(seg, state, jobid);
        job.ToStream(ofs);
        jobid++;
      }
    }
  }

  ofs << "</jobs>" << std::endl;
  ofs.close();
  std::cout << std::endl
            << "... ... In total " << jobid << " jobs" << std::flush;
  return;
}

Job QMMM::createJob(const Segment& seg, const QMState& state, Index jobid) {
  std::string marker = std::to_string(seg.getId()) + ":" + state.ToString();
  std::string tag = seg.getType() + "_" + marker;

  tools::Property Input;
  tools::Property& pInput = Input.add("input", "");
  pInput.add("site_energies", marker);
  tools::Property& regions = pInput.add("regions", "");
  tools::Property& region = regions.add(getFirstRegionName(), "");
  region.add("id", "0");
  if (hasQMRegion()) {
    region.add("state", state.ToString());
  }
  if (use_gs_for_ex_ && (state.Type() == QMStateType::Singlet ||
                         state.Type() == QMStateType::Triplet)) {
    region.add("segments", std::to_string(seg.getId()) + ":n");
  } else {
    region.add("segments", marker);
  }
  Job job(jobid, tag, Input, Job::AVAILABLE);
  return job;
}

void QMMM::ReadJobFile(Topology& top) {

  Index incomplete_jobs = 0;

  Eigen::Matrix<double, Eigen::Dynamic, 5> energies =
      Eigen::Matrix<double, Eigen::Dynamic, 5>::Zero(top.Segments().size(), 5);
  Eigen::Matrix<bool, Eigen::Dynamic, 5> found =
      Eigen::Matrix<bool, Eigen::Dynamic, 5>::Zero(top.Segments().size(), 5);

  tools::Property xml;
  xml.LoadFromXML(jobfile_);
  for (tools::Property* job : xml.Select("jobs.job")) {

    Index jobid = job->get("id").as<Index>();
    if (!job->exists("status")) {
      throw std::runtime_error(
          "Jobfile is malformed. <status> tag missing for job " +
          std::to_string(jobid));
    }
    if (job->get("status").as<std::string>() != "COMPLETE" ||
        !job->exists("output")) {
      incomplete_jobs++;
      continue;
    }

    std::vector<std::string> split =
        tools::Tokenizer(job->get("input.site_energies").as<std::string>(), ":")
            .ToVector();

    Index segid = std::stoi(split[0]);
    if (segid < 0 || segid >= Index(top.Segments().size())) {
      throw std::runtime_error("JobSegment id" + std::to_string(segid) +
                               " is not in topology for job " +
                               std::to_string(jobid));
    }
    QMState state;
    try {
      state.FromString(split[1]);
    } catch (std::runtime_error& e) {
      std::stringstream message;
      message << e.what() << " for job " << jobid;
      throw std::runtime_error(message.str());
    }
    double energy = job->get("output.E_tot").as<double>() * tools::conv::ev2hrt;
    if (found(segid, state.Type().Type()) != 0) {
      throw std::runtime_error("There are two entries in jobfile for segment " +
                               std::to_string(segid) +
                               " state:" + state.ToString());
    }

    energies(segid, state.Type().Type()) = energy;
    found(segid, state.Type().Type()) = true;
  }

  Eigen::Matrix<Index, 1, 5> found_states = found.colwise().count();
  std::cout << std::endl;
  for (Index i = 0; i < 5; i++) {
    if (found_states(i) > 0) {
      QMStateType type(static_cast<QMStateType::statetype>(i));
      std::cout << "Found " << found_states(i) << " states of type "
                << type.ToString() << std::endl;
    }
  }
  if (incomplete_jobs > 0) {
    std::cout << incomplete_jobs << " incomplete jobs found." << std::endl;
  }

  for (Segment& seg : top.Segments()) {
    Index segid = seg.getId();
    for (Index i = 0; i < 4; i++) {
      QMStateType type(static_cast<QMStateType::statetype>(i));
      if (found(segid, i) && found(segid, 4)) {
        double energy = energies(segid, i) - energies(segid, 4);
        seg.setEMpoles(type, energy);
      }
    }
  }
}

}  // namespace xtp
}  // namespace votca
