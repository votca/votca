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
#include <votca/xtp/sqlapplication.h>
#include <votca/xtp/version.h>

namespace votca {
namespace xtp {

SqlApplication::SqlApplication() { Calculatorfactory::RegisterAll(); }

void SqlApplication::Initialize(void) {
  XtpApplication::Initialize();

  Calculatorfactory::RegisterAll();

  namespace propt = boost::program_options;

  AddProgramOptions()("file,f", propt::value<string>(),
                      "  sqlight state file, *.sql");
  AddProgramOptions()("first-frame,i", propt::value<int>()->default_value(1),
                      "  start from this frame");
  AddProgramOptions()("nframes,n", propt::value<int>()->default_value(1),
                      "  number of frames to process");
  AddProgramOptions()("nthreads,t", propt::value<int>()->default_value(1),
                      "  number of threads to create");
  AddProgramOptions()("save,s", propt::value<int>()->default_value(1),
                      "  whether or not to save changes to state file");
}

bool SqlApplication::EvaluateOptions(void) {
  CheckRequired("file", "Please provide the state file");
  return true;
}

void SqlApplication::Run() {

  std::string name = ProgramName();
  if (VersionString() != "") name = name + ", version " + VersionString();
  votca::xtp::HelpTextHeader(name);
  // EVALUATE OPTIONS
  int nThreads = OptionsMap()["nthreads"].as<int>();
  int nframes = OptionsMap()["nframes"].as<int>();
  int fframe = OptionsMap()["first-frame"].as<int>();
  if (fframe-- == 0)
    throw runtime_error(
        "ERROR: First frame is 0, counting "
        "in VOTCA::XTP starts from 1.");
  int save = OptionsMap()["save"].as<int>();

  // STATESAVER & PROGRESS OBSERVER
  string statefile = OptionsMap()["file"].as<string>();
  StateSaverSQLite statsav;
  statsav.Open(_top, statefile);

  // INITIALIZE & RUN CALCULATORS
  cout << "Initializing calculators " << endl;
  BeginEvaluate(nThreads);

  int frameId = -1;
  int framesDone = 0;
  while (statsav.NextFrame() && framesDone < nframes) {
    frameId += 1;
    if (frameId < fframe) continue;
    cout << "Evaluating frame " << _top.getDatabaseId() << endl;
    EvaluateFrame();
    if (save == 1) {
      statsav.WriteFrame();
    } else {
      cout << "Changes have not been written to state file." << endl;
    }
    framesDone += 1;
  }

  if (framesDone == 0)
    cout << "Input requires first frame = " << fframe + 1
         << ", # frames = " << nframes << " => No frames processed.";

  statsav.Close();
  EndEvaluate();
}

void SqlApplication::AddCalculator(ctp::QMCalculator* calculator) {
  _calculators.push_back(calculator);
}

void SqlApplication::BeginEvaluate(int nThreads = 1) {
  for (ctp::QMCalculator* calculator : _calculators) {
    cout << "... " << calculator->Identify() << " ";
    calculator->setnThreads(nThreads);
    calculator->Initialize(&_options);
    cout << endl;
  }
}

bool SqlApplication::EvaluateFrame() {
  for (ctp::QMCalculator* calculator : _calculators) {
    cout << "... " << calculator->Identify() << " " << flush;
    calculator->EvaluateFrame(&_top);
    cout << endl;
  }
  return true;
}

void SqlApplication::EndEvaluate() {
  for (ctp::QMCalculator* calculator : _calculators) {
    calculator->EndEvaluate(&_top);
  }
}

}  // namespace xtp
}  // namespace votca
