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

// Third party includes
#include <boost/format.hpp>

// Local VOTCA includes
#include "votca/xtp/calculatorfactory.h"
#include "votca/xtp/stateapplication.h"
#include "votca/xtp/statesaver.h"
#include "votca/xtp/version.h"

namespace votca {
namespace xtp {

void StateApplication::AddCommandLineOptions() {

  namespace propt = boost::program_options;

  AddProgramOptions()("file,f", propt::value<std::string>(),
                      "  hdf5 state file, *.hdf5");
  AddProgramOptions()("first-frame,i", propt::value<Index>()->default_value(0),
                      "  start from this frame");
  AddProgramOptions()("nframes,n", propt::value<Index>()->default_value(1),
                      "  number of frames to process");
  AddProgramOptions()("save,s", propt::value<bool>()->default_value(true),
                      "  whether or not to save changes to state file");

  AddCommandLineOpt();
}

void StateApplication::EvaluateSpecificOptions() {
  CheckRequired("file", "Please provide the state file");
  CheckRequired("options",
                "Please provide an xml file with calculator options");
  CheckOptions();
}

void StateApplication::execute() {

  options_.LoadFromXML(OptionsMap()["options"].as<std::string>());
  Index nframes = OptionsMap()["nframes"].as<Index>();
  Index fframe = OptionsMap()["first-frame"].as<Index>();
  bool save = OptionsMap()["save"].as<bool>();

  // STATESAVER & PROGRESS OBSERVER
  std::string statefile = OptionsMap()["file"].as<std::string>();
  StateSaver statsav(statefile);
  std::vector<Index> frames = statsav.getFrames();
  if (frames.empty()) {
    throw std::runtime_error("Statefile " + statefile + " not found.");
  }
  // INITIALIZE & RUN CALCULATORS
  std::cout << "Initializing calculator" << std::endl;
  ConfigCalculator();
  std::cout << frames.size() << " frames in statefile, Ids are: ";
  for (Index frame : frames) {
    std::cout << frame << " ";
  }
  std::cout << std::endl;
  if (fframe < Index(frames.size())) {
    std::cout << "Starting at frame " << frames[fframe] << std::endl;
  } else {
    std::cout << "First frame:" << fframe
              << " is larger than number of frames:" << Index(frames.size())
              << std::endl;
    return;
  }

  if ((fframe + nframes) > Index(frames.size())) {
    nframes = Index(frames.size()) - fframe;
  }

  for (Index i = fframe; i < nframes; i++) {
    std::cout << "Evaluating frame " << frames[i] << std::endl;
    Topology top = statsav.ReadFrame(frames[i]);
    EvaluateFrame(top);
    if (save && savetoStateFile()) {
      statsav.WriteFrame(top);
    } else {
      std::cout << "Changes have not been written to state file." << std::endl;
    }
  }
}

}  // namespace xtp
}  // namespace votca
