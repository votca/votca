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
#include <votca/xtp/calculatorfactory.h>
#include <votca/xtp/stateapplication.h>
#include <votca/xtp/version.h>

#include "votca/xtp/statesaver.h"

namespace votca {
namespace xtp {

StateApplication::StateApplication() { Calculatorfactory::RegisterAll(); }

void StateApplication::Initialize(void) {
  XtpApplication::Initialize();

  Calculatorfactory::RegisterAll();

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
}

bool StateApplication::EvaluateOptions(void) {
  CheckRequired("file", "Please provide the state file");
  return true;
}

void StateApplication::Run() {

  std::string name = ProgramName();
  if (VersionString() != "") name = name + ", version " + VersionString();
  xtp::HelpTextHeader(name);

  load_property_from_xml(_options, _op_vm["options"].as<std::string>());
  int nThreads = OptionsMap()["nthreads"].as<int>();
  int nframes = OptionsMap()["nframes"].as<int>();
  int fframe = OptionsMap()["first-frame"].as<int>();
  bool save = OptionsMap()["save"].as<bool>();

  // STATESAVER & PROGRESS OBSERVER
  std::string statefile = OptionsMap()["file"].as<std::string>();
  StateSaver statsav(statefile);
  std::vector<int> frames = statsav.getFrames();
  if (frames.empty()) {
    throw std::runtime_error("Statefile " + statefile + " not found.");
  }
  // INITIALIZE & RUN CALCULATORS
  std::cout << "Initializing calculator" << std::endl;
  BeginEvaluate(nThreads);
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
    std::cout << "Evaluating frame " << frames[i] << std::endl;
    Topology top = statsav.ReadFrame(frames[i]);
    EvaluateFrame(top);
    if (save && _calculator->WriteToStateFile()) {
      statsav.WriteFrame(top);
    } else {
      std::cout << "Changes have not been written to state file." << std::endl;
    }
  }
}

void StateApplication::SetCalculator(QMCalculator* calculator) {
  _calculator = std::unique_ptr<QMCalculator>(calculator);
}

void StateApplication::BeginEvaluate(int nThreads = 1) {
  std::cout << "... " << _calculator->Identify() << " ";
  _calculator->setnThreads(nThreads);
  _calculator->Initialize(_options);
  std::cout << std::endl;
}

bool StateApplication::EvaluateFrame(Topology& top) {
  std::cout << "... " << _calculator->Identify() << " " << std::flush;
  _calculator->EvaluateFrame(top);
  std::cout << std::endl;
  return true;
}

}  // namespace xtp
}  // namespace votca
