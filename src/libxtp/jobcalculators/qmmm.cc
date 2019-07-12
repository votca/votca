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
    jobtop.WriteToHdf5(workdir + "/" + checkpointfilename);

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
  for (const std::unique_ptr<Region>& reg : jobtop) {
    reg->AddResults(regionsresults);
    etot += reg->Etotal();
  }
  std::chrono::time_point<std::chrono::system_clock> end =
      std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;
  jobresult.add("E_tot", std::to_string(etot * tools::conv::hrt2ev));
  jobresult.add("Compute_Time", std::to_string(int(elapsed_time.count())));
  if (!no_top_scf) {
    jobresult.add("Iterations", std::to_string(iteration));
  }
  jres.setOutput(results);
  return jres;
}
void QMMM::WriteJobFile(const Topology& top) {}
void QMMM::ReadJobFile(Topology& top) {}

}  // namespace xtp
};  // namespace votca
