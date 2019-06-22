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
#include <votca/xtp/jobtopology.h>

namespace votca {
namespace xtp {

void QMMM::Initialize(tools::Property& options) {

  std::string key = "options." + Identify();

  _jobfile = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".job_file");

  _print_regions_pdb = options.ifExistsReturnElseReturnDefault(
      key + ".print_regions_pdb", _print_regions_pdb);

  _max_iterations = options.ifExistsReturnElseReturnDefault(
      key + ".max_iterations", _max_iterations);

  if (options.exists(key + ".regions")) {
    _regions_def = options.get(key + ".regions");
  } else {
    throw std::runtime_error("No region definitions found in optionsfile");
  }
}

Job::JobResult QMMM::EvalJob(Topology& top, Job& job, QMThread& Thread) {

  Job::JobResult jres = Job::JobResult();
  Logger& pLog = Thread.getLogger();
  JobTopology jobtop = JobTopology(job, pLog);
  jobtop.BuildRegions(top, _regions_def);
  std::string qmmm_work_dir = "QMMM";
  std::string frame_dir =
      "frame_" + boost::lexical_cast<std::string>(top.getStep());
  if (_print_regions_pdb) {
    std::string pdb_filename =
        "jobtopology_job_" + std::to_string(job.getId()) + ".pdb";
    XTP_LOG_SAVE(logINFO, pLog)
        << "Writing jobtopology to " << pdb_filename << std::flush;
    jobtop.WriteToPdb(pdb_filename);
  }

  int no_static_regions = 0;
  for (std::unique_ptr<Region>& region : jobtop) {
    no_static_regions += region->Converged();
  }
  bool no_top_scf = false;
  if (jobtop.size() - no_static_regions < 2) {
    XTP_LOG_SAVE(logINFO, pLog)
        << "Only " << jobtop.size() - no_static_regions
        << " scf region is used. The remaining regions are static. So no "
           "inter regions scf is required. "
        << std::flush;
    no_top_scf = true;
    _max_iterations = 1;
  }

  for (int iteration = 0; iteration < _max_iterations; iteration++) {

    XTP_LOG_SAVE(logINFO, pLog)
        << " Inter Region SCF Iteration " << iteration + 1 << " of "
        << _max_iterations << std::flush;

    // reverse iterator over regions because the cheapest regions have to be
    // evaluated first
    for (std::vector<std::unique_ptr<Region>>::iterator reg_pointer =
             jobtop.end();
         reg_pointer-- != jobtop.begin();) {
      std::unique_ptr<Region>& region = *reg_pointer;
      XTP_LOG_SAVE(logINFO, pLog)
          << "Running calculations for " << region->identify() << " "
          << region->getId() << std::flush;
      region->Reset();
      region->Evaluate(jobtop.Regions());
    }

    if (!no_top_scf) {
      std::vector<bool> converged_regions;
      for (std::unique_ptr<Region>& region : jobtop) {
        converged_regions.push_back(region->Converged());
      }

      bool all_regions_converged =
          std::all_of(converged_regions.begin(), converged_regions.end(),
                      [](bool i) { return i; });

      if (all_regions_converged) {
        XTP_LOG_SAVE(logINFO, pLog) << "Job converged after " << iteration + 1
                                    << " iterations." << std::flush;
        break;
      }
      if (iteration == _max_iterations - 1) {
        XTP_LOG_SAVE(logINFO, pLog)
            << "Job did not converge after " << iteration + 1
            << " iterations.\n Writing results to jobfile." << std::flush;
      }
    }
  }

  return Job::JobResult();
}
void QMMM::WriteJobFile(Topology& top) {}
void QMMM::ReadJobFile(Topology& top) {}

}  // namespace xtp
};  // namespace votca
