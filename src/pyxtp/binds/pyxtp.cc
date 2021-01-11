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
#include <pybind11/pybind11.h>
#include <vector>

using namespace votca;

/**
 * @brief Construct a new pybind11 module object to invoke votca-xtp*
 *
 */
int call_gwbse(std::string filename) {
  xtp::QMToolFactory::RegisterAll();
  for (const auto& x : xtp::QMTools().getKeys()) {
    std::cout << "key: " << x << "\n";
  }
  // votca::tools::Property prop;
  // prop.LoadFromXML(filename);
  // std::unique_ptr<xtp::QMTool> _tool;
  // _tool = move((xtp::QMTools().Create("dftgwbse")));
  // std::cout << _tool->Identify() << "\n";
  /* std::cout << _tool->Identify(); */
  return 42;
}

PYBIND11_MODULE(xtp_binds, module) {
  module.doc() =
      "VOTCA-XTP is a library which allows you to calculate the electronic "
      "properties of organic materials,"
      "https://votca.github.io";

  module.def("call_calculator", &call_gwbse, R"pbdoc(
        Invoke a Votca XTP calculator

        Parameters
        ----------
        name
          Calculator's name

  )pbdoc");
}
