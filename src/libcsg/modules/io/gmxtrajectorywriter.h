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

#ifndef _VOTCA_CSG_GMXTRAJECTORYWRITER_H
#define _VOTCA_CSG_GMXTRAJECTORYWRITER_H

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include <votca/csg/topology.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/unitconverter.h>

struct t_trxstatus;
namespace votca {
namespace csg {

class GMXTrajectoryWriter : public TrajectoryWriter {
 public:
  const tools::DistanceUnit distance_unit = tools::DistanceUnit::nanometers;
  const tools::MassUnit mass_unit = tools::MassUnit::atomic_mass_units;
  const tools::TimeUnit time_unit = tools::TimeUnit::picoseconds;
  const tools::ChargeUnit charge_unit = tools::ChargeUnit::e;
  const tools::MolarEnergyUnit energy_unit =
      tools::EnergyUnit::kilojoules_per_mole;
  const tools::VelocityUnit velocity_unit =
      tools::VelocityUnit::nanometers_per_picosecond;
  const tools::MolarForceUnit force_unit =
      tools::ForceUnit::kilojoules_per_mole_nanometer;

  GMXTrajectoryWriter() = default;

  void Open(std::string file, bool = false) override;
  void Close() override;
  void Write(Topology *conf) override;

 private:
  t_trxstatus *_file;
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_GMXTRAJECTORYWRITER_H */
