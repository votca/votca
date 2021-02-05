/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "votca/xtp/statesaver.h"
#include "xtp_bind_calculators.h"
#include <iostream>
#include <pybind11/embed.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace votca;

namespace pyxtp {

void call_calculator(const std::string& calculatorName,
                     const std::map<std::string, std::string>& dict) {

  // Retrieve arguments from Python dictionary
  const std::string& xmlFile = dict.at("xmlFile");
  const std::string& stateFile = dict.at("stateFile");
  Index nThreads = std::stol(dict.at("nThreads"));
  Index firstFrame = std::stol(dict.at("firstFrame"));
  Index nFrames = std::stol(dict.at("nFrames"));
  bool save = boost::lexical_cast<bool>("0");

  // Load properties
  votca::tools::Property prop;
  prop.LoadFromXML(xmlFile);
  // Call calculator
  pyxtp::XTPCalculators calc;
  calc.Initialize(calculatorName, nThreads, prop);
  calc.Run(stateFile, nFrames, firstFrame, save);
}

void XTPCalculators::Initialize(const std::string& name, Index nThreads,
                                votca::tools::Property prop) {
  xtp::Calculatorfactory factory;
  _calculator = factory.Create(name);
  _calculator->setnThreads(nThreads);
  _calculator->Initialize(prop);
  std::cout << "Calculator has been Initialized\n";
}

void XTPCalculators::Run(const std::string& stateFile, Index nFrames,
                         Index firstFrame, bool save) {
  // STATESAVER & PROGRESS OBSERVER
  xtp::StateSaver statsav(stateFile);
  std::vector<Index> frames = statsav.getFrames();
  if (frames.empty()) {
    throw std::runtime_error("stateFile " + stateFile + " not found.");
  }
  // INITIALIZE & RUN CALCULATORS
  std::cout << "Initializing calculator\n";
  std::cout << frames.size() << " frames in stateFile, Ids are: ";
  for (Index frame : frames) {
    std::cout << frame << " ";
  }
  std::cout << "\n";
  if (firstFrame < Index(frames.size())) {
    std::cout << "Starting at frame " << frames[firstFrame] << std::endl;
  } else {
    std::cout << "First frame:" << firstFrame
              << " is larger than number of frames:" << Index(frames.size())
              << std::endl;
    return;
  }

  if ((firstFrame + nFrames) > Index(frames.size())) {
    nFrames = Index(frames.size()) - firstFrame;
  }

  for (Index i = firstFrame; i < nFrames; i++) {
    std::cout << "Evaluating frame " << frames[i] << std::endl;
    xtp::Topology top = statsav.ReadFrame(frames[i]);
    _calculator->EvaluateFrame(top);
    if (save && _calculator->WriteToStateFile()) {
      statsav.WriteFrame(top);
    } else {
      std::cout << "Changes have not been written to state file." << std::endl;
    }
  }
}

}  // namespace pyxtp
