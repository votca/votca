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

#include <boost/format.hpp>
#include <votca/xtp/jobapplication.h>
#include <votca/xtp/jobcalculatorfactory.h>
#include <votca/xtp/version.h>

namespace votca {
namespace xtp {

JobApplication::JobApplication() { JobCalculatorfactory::RegisterAll(); }

void JobApplication::Initialize(void) {
  XtpApplication::Initialize();

  JobCalculatorfactory::RegisterAll();

  namespace propt = boost::program_options;

  AddProgramOptions()("file,f", propt::value<std::string>(),
                      "  sqlite state file, *.sql");
  AddProgramOptions()("first-frame,i", propt::value<int>()->default_value(1),
                      "  start from this frame");
  AddProgramOptions()("nframes,n", propt::value<int>()->default_value(1),
                      "  number of frames to process");
  AddProgramOptions()("nthreads,t", propt::value<int>()->default_value(1),
                      "  number of threads to create");
  AddProgramOptions()("save,s", propt::value<int>()->default_value(1),
                      "  whether or not to save changes to state file");
  AddProgramOptions()("restart,r",
                      propt::value<std::string>()->default_value(""),
                      "  restart pattern: 'host(pc1:234) stat(FAILED)'");
  AddProgramOptions()("cache,c", propt::value<int>()->default_value(8),
                      "  assigns jobs in blocks of this size");
  AddProgramOptions()("jobs,j",
                      propt::value<std::string>()->default_value("run"),
                      "  task(s) to perform: input, run, import");
  AddProgramOptions()("maxjobs,m", propt::value<int>()->default_value(-1),
                      "  maximum number of jobs to process (-1 = inf)");
}

bool JobApplication::EvaluateOptions(void) {
  CheckRequired("options",
                "Please provide an xml file with calculator options");
  CheckRequired("file", "Please provide the state file");

  std::string jobstr = _op_vm["jobs"].as<std::string>();
  _generate_input = jobstr.find("write") != std::string::npos;
  _run = jobstr.find("run") != std::string::npos;
  _import = jobstr.find("read") != std::string::npos;

  return true;
}

void JobApplication::Run() {
  std::string name = ProgramName();
  if (VersionString() != "") name = name + ", version " + VersionString();
  HelpTextHeader(name);

  load_property_from_xml(_options, _op_vm["options"].as<std::string>());

  // EVALUATE OPTIONS
  int nThreads = OptionsMap()["nthreads"].as<int>();
  int nframes = OptionsMap()["nframes"].as<int>();
  int fframe = OptionsMap()["first-frame"].as<int>();
  if (fframe-- == 0)
    throw std::runtime_error(
        "ERROR: First frame is 0, counting "
        "in VOTCA::XTP starts from 1.");
  int save = OptionsMap()["save"].as<int>();

  // STATESAVER & PROGRESS OBSERVER
  std::string statefile = OptionsMap()["file"].as<std::string>();
  StateSaverSQLite statsav;
  statsav.Open(statefile);

  ProgObserver<std::vector<Job*>, Job*, Job::JobResult> progObs =
      ProgObserver<std::vector<Job*>, Job*, Job::JobResult>();
  progObs.InitCmdLineOpts(OptionsMap());

  // INITIALIZE & RUN CALCULATORS
  std::cout << "Initializing calculators " << std::endl;
  BeginEvaluate(nThreads, &progObs);

  int frameId = -1;
  int framesDone = 0;
  while (statsav.NextFrame() && framesDone < nframes) {
    frameId += 1;
    if (frameId < fframe) continue;
    std::cout << "Evaluating frame " << _top.getStep() << std::endl;
    EvaluateFrame();
    if (save == 1) {
      statsav.WriteFrame(_top);
    } else {
      std::cout << "Changes have not been written to state file." << std::endl;
    }
    framesDone += 1;
  }

  if (framesDone == 0)
    std::cout << "Input requires first frame = " << fframe + 1
              << ", # frames = " << nframes << " => No frames processed.";

  statsav.Close();
  EndEvaluate();
}

void JobApplication::AddCalculator(JobCalculator* calculator) {
  _calculators.push_back(std::unique_ptr<JobCalculator>(calculator));
}

void JobApplication::BeginEvaluate(
    int nThreads = 1,
    ProgObserver<std::vector<Job*>, Job*, Job::JobResult>* obs = NULL) {

  for (std::unique_ptr<JobCalculator>& calculator : _calculators) {
    std::cout << "... " << calculator->Identify() << " ";
    calculator->setnThreads(nThreads);
    calculator->setProgObserver(obs);
    calculator->Initialize(_options);
    std::cout << std::endl;
  }
}

bool JobApplication::EvaluateFrame() {
  for (std::unique_ptr<JobCalculator>& calculator : _calculators) {
    std::cout << "... " << calculator->Identify() << " " << std::flush;
    if (_generate_input) calculator->WriteJobFile(_top);
    if (_run) calculator->EvaluateFrame(_top);
    if (_import) calculator->ReadJobFile(_top);
    std::cout << std::endl;
  }
  return true;
}

void JobApplication::EndEvaluate() {
  for (std::unique_ptr<JobCalculator>& calculator : _calculators) {
    calculator->EndEvaluate(_top);
  }
}

}  // namespace xtp
}  // namespace votca
