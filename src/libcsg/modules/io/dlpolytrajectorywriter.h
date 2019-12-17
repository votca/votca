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

#ifndef _DLPOLYTRAJECTORYWRITER_H
#define _DLPOLYTRAJECTORYWRITER_H

#include <votca/csg/topology.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/unitconverter.h>

namespace votca {
namespace csg {

/**
    \brief class for writing dlpoly trajectory and configuration files

    This class encapsulates the dlpoly trajectory and configuration writing
   function

*/

class DLPOLYTrajectoryWriter : public TrajectoryWriter {
 public:
  const tools::DistanceUnit distance_unit = tools::DistanceUnit::angstroms;
  const tools::MassUnit mass_unit = tools::MassUnit::atomic_mass_units;
  const tools::TimeUnit time_unit = tools::TimeUnit::picoseconds;
  const tools::ChargeUnit charge_unit = tools::ChargeUnit::e;
  const tools::MolarEnergyUnit energy_unit = tools::EnergyUnit::joules_per_mole;
  const tools::VelocityUnit velocity_unit =
      tools::VelocityUnit::angstroms_per_picosecond;
  const tools::MolarForceUnit force_unit =
      tools::ForceUnit::kilojoules_per_mole_angstrom;

  // open transformed trajectory file
  void Open(std::string file, bool bAppend = false) override;
  // close transformed trajectory file
  void Close() override;
  // write a frame into transformed trajectory file
  void Write(Topology *conf) override;

  /// set/get the created configuration or trajectory file name:
  /// <name>.dlpc or <name>.dlph (convention: ".dlpc"="CONFIG_CGV",
  /// ".dlph"="HISTORY_CGV")
  void setFname(std::string name) {
    _fname = name;
    return;
  }
  std::string getFname() { return _fname; }

  /// set/check the flag for the created file as configuration, i.e. not
  /// trajectory format
  void setIsConfig(bool isConf) {
    _isConfig = isConf;
    return;
  }
  bool getIsConfig() { return _isConfig; }

 private:
  std::ofstream _fl;
  std::string _fname;
  bool _isConfig;
};

}  // namespace csg
}  // namespace votca

#endif /* _DLPOLYTRAJECTORYWRITER_H */
