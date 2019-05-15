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

  if (options.exists(key + ".regions")) {
    _regions_def = options.get(key + ".regions");
  } else {
    throw std::runtime_error("No region definitions found in optionsfile");
  }

  if (options.exists(key + ".interactors")) {
    _interactor_def = options.get(key + ".interactors");
  } else {
    throw std::runtime_error("No interactor definitions found in optionsfile");
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
  return Job::JobResult();
}
void QMMM::WriteJobFile(Topology& top) {}
void QMMM::ReadJobFile(Topology& top) {}

}  // namespace xtp
};  // namespace votca
