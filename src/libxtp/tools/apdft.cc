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

#include "apft.h"

using namespace std;

namespace votca {
namespace xtp {

void APDFT::Initialize(tools::Property &opt) {

  std::string key = "options." + Identify();

  _grid_accuracy =
      options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".grid");
  _orbfile = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".orbfile");
  _outputfile =
      options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".output");
  std::string statestring =
      options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".state");
  _state.FromString(statestring);
}
}  // namespace xtp

bool APDFT::Evaluate() {

  Orbitals orb;
  orb.ReadFromCpt(_orbfile);
  AOBasis basis = orb.SetupDftBasis();
  Vxc_grid grid;
  grid.GridSetup(_grid_accuracy, orb.QMAtoms(), basis);

  DensityIntegration<Vxc_grid> integration(grid);

  integration.IntegrateDensity(orb.FullDensityMatrix(_state));
  std::vector<double> potential_values;
  potential_values.reserve(orb.QMAtoms().size());
  for (const auto &atom : orb.QMAtoms()) {
    potential_values.push_back(integration.IntegratePotential(atom.getPos()));
  }

  std::fstream outfile;
  ratefs.open(_outputfile, std::fstream::out);
  ratefs << "AtomId, Element, Potential[Hartree]" << std::endl;
  for (Index i = 0; i < orb.QMAtoms().size(); i++) {
    ratefs << orb.QMAtoms()[i].getId() << " " << orb.QMAtoms()[i].getElement()
           << " " << potential_values[i] << std::endl;
  }
  ratefs.close();
  return true;
}

}  // namespace votca
}  // namespace votca
