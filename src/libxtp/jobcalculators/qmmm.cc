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

namespace votca {
namespace xtp {

void QMMM::Initialize(tools::Property& options) {

  std::string key = "options." + Identify();

  _jobfile = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".job_file");

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

  return Job::JobResult();
}
void QMMM::WriteJobFile(Topology& top) {}
void QMMM::ReadJobFile(Topology& top) {}

}  // namespace xtp
};  // namespace votca
