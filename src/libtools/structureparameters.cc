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
#include "../../include/votca/tools/structureparameters.h"
#include "../../include/votca/tools/unitconverter.h"

namespace votca {
namespace tools {

static UnitConverter converter_;

void StructureParameters::convertParameterIfPossible_(
    const StructureParameter parameter) {
  if (parameter == StructureParameter::XTP_Position) {
    if (parameters_.count(StructureParameter::CSG_Position)) {
      parameters_[StructureParameter::XTP_Position] =
          converter_.convert(DistanceUnit::nanometers,
                             DistanceUnit::angstroms) *
          boost::any_cast<Eigen::VectorXd>(
              parameters_[StructureParameter::CSG_Position]);
    }
  } else if (parameter == CSG_Position) {
    if (parameters_.count(StructureParameter::XTP_Position)) {
      parameters_[StructureParameter::CSG_Position] =
          converter_.convert(DistanceUnit::nanometers,
                             DistanceUnit::angstroms) *
          boost::any_cast<Eigen::VectorXd>(
              parameters_[StructureParameter::XTP_Position]);
    }
  }
}

}  // namespace tools
}  // namespace votca
