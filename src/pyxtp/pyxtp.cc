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

// #include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


/**
 * @brief Construct a new pybind11 module object to invoke votca-xtp*
 * 
 */
int call_calculator(int i, int j) {
    return i + j;
}


PYBIND11_MODULE(pyxtp, module) {
  module.doc() =
      "VOTCA-XTP is a library which allows you to calculate the electronic "
      "properties of organic materials,"
      "https://votca.github.io";

  module.def("call_calculator", &call_calculator, R"pbdoc(
        Invoke a Votca XTP calculator

        Parameters
        ----------
        name
          Calculator's name

  )pbdoc");
}
