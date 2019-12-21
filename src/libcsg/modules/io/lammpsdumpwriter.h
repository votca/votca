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

#ifndef __VOTCA_CSG_LAMMPSDUMPWRITER_H
#define __VOTCA_CSG_LAMMPSDUMPWRITER_H

#include "../../../../include/votca/csg/topology.h"
#include "../../../../include/votca/csg/trajectorywriter.h"
#include <stdio.h>
#include <votca/tools/unitconverter.h>

namespace votca {
namespace csg {

class LAMMPSDumpWriter : public TrajectoryWriter {
 public:
  /// Assuming units are using 'units real' lammps command
  const tools::DistanceUnit distance_unit = tools::DistanceUnit::angstroms;
  const tools::TimeUnit time_unit = tools::TimeUnit::femtoseconds;
  const tools::MassUnit mass_unit = tools::MassUnit::grams_per_mole;
  const tools::MolarEnergyUnit energy_unit =
      tools::MolarEnergyUnit::kilocalories_per_mole;
  const tools::ChargeUnit charge_unit = tools::ChargeUnit::e;
  const tools::MolarForceUnit force_unit =
      tools::MolarForceUnit::kilocalories_per_mole_angstrom;
  const tools::VelocityUnit velocity_unit =
      tools::VelocityUnit::angstroms_per_femtosecond;

  void Open(std::string file, bool bAppend = false) override;
  void Close() override;

  void Write(Topology *conf) override;

 private:
  FILE *_out;
};

}  // namespace csg
}  // namespace votca

#endif  // __VOTCA_CSG_LAMMPSDUMPWRITER_H
