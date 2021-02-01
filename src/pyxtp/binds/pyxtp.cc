/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#include "pyxtp.hpp"
#include <iostream>
#include <memory>
#include <pybind11/stl.h>
#include <string>
#include <vector>

using namespace votca;

namespace pyxtp {

/**
 * @brief Construct a new pybind11 module object to invoke votca-xtp*
 *
 */
int call_calculator(const std::string& name, int n_threads,
                    std::string xml_file) {
  // Load properties
  votca::tools::Property prop;
  prop.LoadFromXML(xml_file);
  // Call calculator
  pyxtp::PyXTP pyxtp;
  pyxtp.Initialize(name, n_threads, prop);
  return 42;
}

void PyXTP::Initialize(const std::string& name, int n_threads,
                       votca::tools::Property prop) {
  xtp::Calculatorfactory inst;
  _calculator = inst.Create(name);
  _calculator->setnThreads(n_threads);
  _calculator->Initialize(prop);
  std::cout << "The calculator has been initialize!\n";
}

/**
 * @brief Construct a Votca tools property using the provided dictionary
 *
 */

}  // namespace pyxtp

PYBIND11_MODULE(xtp_binds, module) {
  module.doc() =
      "VOTCA-XTP is a library which allows you to calculate the electronic "
      "properties of organic materials,"
      "https://votca.github.io";

  module.def("call_calculator", &pyxtp::call_calculator,
             R"pbdoc(
        Invoke a Votca XTP calculator

        Parameters
        ----------
        name
          Calculator's name

  )pbdoc");
}
