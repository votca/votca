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

#include "xtp_bind_calculators.h"
#include <iostream>
#include <pybind11/embed.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace votca;

namespace pyxtp {

/**
 * @brief Construct a new pybind11 module object to invoke votca-xtp*
 *
 */
int call_calculator(const std::string& name,
                    const std::map<std::string, std::string>& dict) {
  // Load properties
  votca::tools::Property prop;
  const std::string& xml_file = dict.at("xml_file");
  prop.LoadFromXML(xml_file);
  // Call calculator
  pyxtp::XTPCalculators calc;
  calc.Initialize(name, 1, prop);
  return 42;
}

void XTPCalculators::Initialize(const std::string& name, int n_threads,
                                votca::tools::Property prop) {
  xtp::Calculatorfactory factory;
  _calculator = factory.Create(name);
  _calculator->setnThreads(n_threads);
  _calculator->Initialize(prop);
  std::cout << "The calculator has been initialize!\n";
}

}  // namespace pyxtp
