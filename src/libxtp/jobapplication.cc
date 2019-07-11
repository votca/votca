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
                      "  hdf5 state file, *.hdf5");
  AddProgramOptions()("first-frame,i", propt::value<int>()->default_value(0),
                      "  start from this frame");
  AddProgramOptions()("nframes,n", propt::value<int>()->default_value(1),
                      "  number of frames to process");
  AddProgramOptions()("nthreads,t", propt::value<int>()->default_value(1),
                      "  number of threads to create");
  AddProgramOptions()("save,s", propt::value<bool>()->default_value(true),
                      "  whether or not to save changes to state file");
  AddProgramOptions()("restart,r",
                      propt::value<std::string>()->default_value(""),
                      "  restart pattern: 'host(pc1:234) stat(FAILED)'");
  AddProgramOptions()("cache,c", propt::value<int>()->default_value(8),
                      "  assigns jobs in blocks of this size");
  AddProgramOptions()("jobs,j",
                      propt::value<std::string>()->default_value("run"),
                      "  task(s) to perform: write, run, read");
  AddProgramOptions()("maxjobs,m", propt::value<int>()->default_value(-1),
                      "  maximum number of jobs to process (-1 = inf)");
}

bool JobApplication::EvaluateOptions(void) {
  CheckRequired("options",
                "Please provide an xml file with calculator options");
  CheckRequired("file", "Please provide the state file");

  std::string jobstr = _op_vm["jobs"].as<std::string>();
  _generate_input = (jobstr == "write");
  _run = (jobstr == "run");
  _import = (jobstr == "read");

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
  bool save = OptionsMap()["save"].as<bool>();

  // STATESAVER & PROGRESS OBSERVER
  std::string statefile = OptionsMap()["file"].as<std::string>();
  ProgObserver<std::vector<Job>> progObs = ProgObserver<std::vector<Job>>();
  progObs.InitCmdLineOpts(OptionsMap());

  // INITIALIZE & RUN CALCULATORS
  std::cout << "Initializing calculator " << std::endl;
  BeginEvaluate(nThreads, progObs);

  StateSaver statsav(statefile);

  std::vector<int> frames = statsav.getFrames();

  std::cout << frames.size() << " frames in statefile, Ids are: ";
  for (int frame : frames) {
    std::cout << frame << " ";
  }
  std::cout << std::endl;
  if (fframe < int(frames.size())) {
    std::cout << "Starting at frame " << frames[fframe] << std::endl;
  } else {
    std::cout << "First frame:" << fframe
              << " is larger than number of frames:" << int(frames.size())
              << std::endl;
    return;
  }

  if ((fframe + nframes) > int(frames.size())) {
    nframes = frames.size() - fframe;
  }

  for (int i = fframe; i < nframes; i++) {
    std::cout << "Evaluating frame " << i << std::endl;
    Topology top = statsav.ReadFrame(i);
    EvaluateFrame(top);
    if (save && _import) {
      statsav.WriteFrame(top);
    } else {
      std::cout << "Changes have not been written to state file." << std::endl;
    }
  }
}

void JobApplication::SetCalculator(JobCalculator* calculator) {
  _calculator = std::unique_ptr<JobCalculator>(calculator);
}

void JobApplication::BeginEvaluate(int nThreads,
                                   ProgObserver<std::vector<Job>>& obs) {

  std::cout << "... " << _calculator->Identify() << " ";
  _calculator->setnThreads(nThreads);
  _calculator->setProgObserver(&obs);
  _calculator->Initialize(_options);
  std::cout << std::endl;
}

bool JobApplication::EvaluateFrame(Topology& top) {
  std::cout << "... " << _calculator->Identify() << " " << std::flush;
  if (_generate_input) {
    _calculator->WriteJobFile(top);
  } else if (_run) {
    _calculator->EvaluateFrame(top);
  } else if (_import) {
    _calculator->ReadJobFile(top);
  } else {
    ;
  }
  std::cout << std::endl;
  return true;
}

}  // namespace xtp
}  // namespace votca
