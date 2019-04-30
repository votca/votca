/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#ifndef VOTCA_TOOLS_PARAMETERS_H
#define VOTCA_TOOLS_PARAMETERS_H

#include <boost/any.hpp>
#include <cassert>
#include <unordered_map>

namespace votca {
namespace tools {

/**
    \brief Provides a means to reduce cross dependency of header files
 *
 * This templace class allows the conversion of one type of object to another
 * this is done by breaking the said object up into the standard types.
 *
 * For the type converter to work each object it is converting between must
 * have defined at least one write function and or one read function.
 *
 * This tempalte object only allows converting between doubles,strings,ints.
 */
class Parameters {

 public:
  enum Parameter {
    Mass,
    Position,
    MoleculeId,
    ResidueId,
    Charge,
    Element,
    Symmetry,
    ResidueType,
    BeadId,
    BeadType,
    MoleculeType
  };

  void set(const Parameter parameter, boost::any value);

  template <class T>
  T get(const Parameter parameter) const;

 private:
  std::unordered_map<Parameter, boost::any> parameters;
};

void Parameters::set(const Parameter parameter, boost::any value) {
  parameters[parameter] = value;
}

template <class T>
T Parameters::get(const Parameter parameter) const {
  assert(parameters.count(parameter) &&
         "Parameter is not stored in Parameters class");
  return parameters.at(parameter);
}

}  // namespace tools
}  // namespace votca
#endif  // # VOTCA_TOOL_PARAMETERS_H
