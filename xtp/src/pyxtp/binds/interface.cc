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
#include "xtp_bind_tools.h"
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

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
        threads
          Number of threads to perform the computation
        xml_file
          Input file specification
  )pbdoc");
  module.def("call_tool", &pyxtp::call_tool,
             R"pbdoc(
        Invoke a Votca XTP Tool

        Parameters
        ----------
        name
          Calculator's name
        threads
          Number of threads to perform the computation
        xml_file
          Input file specification
  )pbdoc");
}
