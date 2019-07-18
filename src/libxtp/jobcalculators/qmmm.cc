/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include "qmmm.h"
#include <boost/filesystem.hpp>
#include <chrono>
#include <votca/xtp/jobtopology.h>

namespace votca {
namespace xtp {

void QMMM::Initialize(tools::Property& options) {

  std::string key = "options." + Identify();
  ParseCommonOptions(options);
  _jobfile = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".job_file");

  _print_regions_pdb = options.ifExistsReturnElseReturnDefault(
      key + ".print_regions_pdb", _print_regions_pdb);

  _max_iterations = options.ifExistsReturnElseReturnDefault(
      key + ".max_iterations", _max_iterations);

  if (options.exists(key + ".regions")) {
    _regions_def = options.get(key + ".regions");
    _regions_def.add("mapfile", _mapfile);
  } else {
    throw std::runtime_error("No region definitions found in optionsfile");
  }

  if (options.exists(key + ".write_parse")) {
    _write_parse = true;

    std::string states = options.ifExistsReturnElseReturnDefault<std::string>(
        key + ".write_parse.states", "e h");
    tools::Tokenizer tok(states, " ,;\n\t");
    std::vector<std::string> statestrings = tok.ToVector();
    _states.reserve(statestrings.size());
    for (std::string s : statestrings) {
      _states.push_back(QMStateType(s));
    }
    bool groundstate_found = false;
    for (const QMStateType& state : _states) {
      if (state.Type() == QMStateType::Gstate) groundstate_found = true;
    }
    if (!groundstate_found) {
      _states.push_back(QMStateType("n"));
    }
  }
}

Job::JobResult QMMM::EvalJob(const Topology& top, Job& job, QMThread& Thread) {
  std::chrono::time_point<std::chrono::system_clock> start =
      std::chrono::system_clock::now();

  std::string qmmm_work_dir = "QMMM";
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
  jobtop.BuildRegions(top, _regions_def);

  if (_print_regions_pdb) {
    std::string pdb_filename = "regions.pdb";
    XTP_LOG_SAVE(logINFO, pLog) << TimeStamp() << " Writing jobtopology to "
                                << (workdir + "/" + pdb_filename) << std::flush;
    jobtop.WriteToPdb(workdir + "/" + pdb_filename);
  }

  int no_static_regions = 0;
  for (std::unique_ptr<Region>& region : jobtop) {
    no_static_regions += region->Converged();
  }
  bool no_top_scf = false;
  if (jobtop.size() - no_static_regions < 2) {
    XTP_LOG_SAVE(logINFO, pLog)
        << TimeStamp() << " Only " << jobtop.size() - no_static_regions
        << " scf region is used. The remaining regions are static. So no "
           "inter regions scf is required. "
        << std::flush;
    no_top_scf = true;
    _max_iterations = 1;
  }
  int iteration = 0;
  for (; iteration < _max_iterations; iteration++) {

    XTP_LOG_SAVE(logINFO, pLog)
        << TimeStamp() << " --Inter Region SCF Iteration " << iteration + 1
        << " of " << _max_iterations << std::flush;

    for (std::unique_ptr<Region>& region : jobtop) {
      XTP_LOG_SAVE(logINFO, pLog)
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
    XTP_LOG_SAVE(logINFO, pLog) << TimeStamp() << " Writing checkpoint to "
                                << checkpointfilename << std::flush;
    // jobtop.WriteToHdf5(workdir + "/" + checkpointfilename);

    if (!no_top_scf) {
      std::vector<bool> converged_regions;
      for (std::unique_ptr<Region>& region : jobtop) {
        converged_regions.push_back(region->Converged());
      }

      double etot = 0.0;
      for (const std::unique_ptr<Region>& reg : jobtop) {
        etot += reg->Etotal();
      }
      XTP_LOG_SAVE(logINFO, pLog)
          << TimeStamp() << " --Total Energy all regions " << etot
          << std::flush;

      bool all_regions_converged =
          std::all_of(converged_regions.begin(), converged_regions.end(),
                      [](bool i) { return i; });

      if (all_regions_converged) {
        XTP_LOG_SAVE(logINFO, pLog)
            << TimeStamp() << " Job converged after " << iteration + 1
            << " iterations." << std::flush;
        jres.setStatus(Job::JobStatus::COMPLETE);
        break;
      }
      if (iteration == _max_iterations - 1) {
        XTP_LOG_SAVE(logINFO, pLog)
            << TimeStamp() << " Job did not converge after " << iteration + 1
            << " iterations.\n Writing results to jobfile." << std::flush;
        jres.setStatus(Job::JobStatus::FAILED);
        jres.setError("Inter Region SCF did not converge in " +
                      std::to_string(_max_iterations) + " iterations.");
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
  jobresult.add("Compute_Time", std::to_string(int(elapsed_time.count())));
  jobresult.add("Total_Charge", std::to_string(charge));
  if (!no_top_scf) {
    jobresult.add("Iterations", std::to_string(iteration + 1));
  }
  jres.setOutput(results);
  return jres;
}
void QMMM::WriteJobFile(const Topology& top) {

  if (!_write_parse) {
    throw std::runtime_error(
        "Cannot write jobfile, please add <write_parse><states>e "
        "s1</states></write_parse> to your options.");
  }

  std::cout << std::endl
            << "... ... Writing job file " << _jobfile << std::flush;

  std::ofstream ofs;
  ofs.open(_jobfile, std::ofstream::out);
  if (!ofs.is_open())
    throw std::runtime_error("\nERROR: bad file handle: " + _jobfile);

  ofs << "<jobs>" << std::endl;
  int jobid = 0;
  for (const Segment& seg : top.Segments()) {
    for (const QMStateType& state : _states) {

      std::string marker = std::to_string(seg.getId()) + ":" + state.ToString();
      std::string tag = seg.getName() + "_" + marker;

      tools::Property Input;
      tools::Property& pInput = Input.add("input", "");
      pInput.add("site_energies", marker);
      tools::Property& regions = pInput.add("regions", "");
      tools::Property& region = regions.add("region", "");
      region.add("id", "0");
      region.add("segments", marker);
      Job job(jobid, tag, Input, Job::AVAILABLE);
      job.ToStream(ofs);
      jobid++;
    }
  }

  ofs << "</jobs>" << std::endl;
  ofs.close();
  std::cout << std::endl
            << "... ... In total " << jobid + 1 << " jobs" << std::flush;
  return;
}
void QMMM::ReadJobFile(Topology& top) {

  if (!_write_parse) {
    throw std::runtime_error(
        "Cannot read jobfile, please add <write_parse><states>n e "
        "h</states></write_parse> to your options.");
  }

  int incomplete_jobs = 0;

  Eigen::Matrix<double, Eigen::Dynamic, 5> energies =
      Eigen::Matrix<double, Eigen::Dynamic, 5>::Zero(top.Segments().size(), 5);
  Eigen::Matrix<int, Eigen::Dynamic, 5> found =
      Eigen::Matrix<int, Eigen::Dynamic, 5>::Zero(top.Segments().size(), 5);

  tools::Property xml;
  load_property_from_xml(xml, _jobfile);
  std::vector<tools::Property*> jobProps = xml.Select("jobs.job");
  for (tools::Property* job : jobProps) {

    int jobid = job->get("id").as<int>();
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

    std::string marker = job->get("input.site_energies").as<std::string>();
    tools::Tokenizer tok(marker, ":");
    std::vector<std::string> split = tok.ToVector();
    int segid = std::stoi(split[0]);
    if (segid < 0 || segid >= int(top.Segments().size())) {
      throw std::runtime_error("JobSegment id" + std::to_string(segid) +
                               " is not in topology for job " +
                               std::to_string(jobid));
    }
    QMStateType state;
    try {
      state.FromString(split[1]);
    } catch (std::runtime_error& e) {
      std::stringstream message;
      message << e.what() << " for job " << jobid;
      throw std::runtime_error(message.str());
    }

    double energy = job->get("output.E_tot").as<double>();
    if (found(segid, state.Type()) != 0) {
      throw std::runtime_error("There are two entries in jobfile for segment " +
                               std::to_string(segid) +
                               " state:" + state.ToString());
    }

    energies(segid, state.Type()) = energy;
    found(segid, state.Type()) = 1;
  }

  Eigen::Matrix<int, 1, 5> found_states = found.colwise().sum();
  for (int i = 0; i < 5; i++) {
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
    int segid = seg.getId();
    for (int i = 0; i < 4; i++) {
      QMStateType type(static_cast<QMStateType::statetype>(i));
      if (found(segid, i) && found(segid, 4)) {
        double energy = energies(segid, i) - energies(segid, 4);
        seg.setEMpoles(type, energy);
      }
    }
  }
}

}  // namespace xtp
};  // namespace votca
