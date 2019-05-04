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
#pragma once
#ifndef VOTCA_TOOLS_UNITCONVERTER_H
#define VOTCA_TOOLS_UNITCONVERTER_H

#include <iostream>
#include <map>
namespace votca {
namespace tools {

enum DistanceUnit { meters, centimeters, nanometers, angstroms, bohr };

enum MassUnit {
  attograms,
  picograms,
  femtograms,
  atomic_mass_units,
  kilograms,
  grams
};

enum TimeUnit { seconds, microseconds, nanoseconds, femtoseconds, picoseconds };

enum EnergyUnit { electron_volts, kilocalories, hartrees, joules };

/**
 * @brief Class converts between different unit types
 */
class UnitConverter {
 private:
  /// All distances with respect to Ang
  const std::map<DistanceUnit, double> distance_map_ = {
      {DistanceUnit::meters, 1E10},
      {DistanceUnit::centimeters, 1E8},
      {DistanceUnit::nanometers, 10.0},
      {DistanceUnit::angstroms, 1.0},
      {DistanceUnit::bohr, 0.529177}};

  // All times with respect to pico seconds
  const std::map<TimeUnit, double> time_map_ = {
      {TimeUnit::seconds, 1E12},
      {TimeUnit::microseconds, 1E6},
      {TimeUnit::nanoseconds, 1000},
      {TimeUnit::picoseconds, 1.0},
      {TimeUnit::femtoseconds, 0.001}};

  // All masses in terms of atomic mass units
  const std::map<MassUnit, double> mass_map_ = {
      {MassUnit::kilograms, 6.02214E26}, {MassUnit::grams, 6.02214E23},
      {MassUnit::picograms, 6.02214E14}, {MassUnit::femtograms, 6.02214E11},
      {MassUnit::attograms, 602214},     {MassUnit::atomic_mass_units, 1.0}};

  // All energies in terms of electron volts
  const std::map<EnergyUnit, double> energy_map_ = {
      {EnergyUnit::kilocalories, 2.613195131836172E22},
      {EnergyUnit::joules, 6.242E18},
      {EnergyUnit::hartrees, 27.2114},
      {EnergyUnit::electron_volts, 1.0}};

 public:
  double convert(const DistanceUnit& from, const DistanceUnit& to,
                 const double& from_value) const noexcept {
    return distance_map_.at(from) / distance_map_.at(to) * from_value;
  }
  double convert(const TimeUnit& from, const TimeUnit& to,
                 const double& from_value) const noexcept {
    return time_map_.at(from) / time_map_.at(to) * from_value;
  }
  double convert(const MassUnit& from, const MassUnit& to,
                 const double& from_value) const noexcept {
    return mass_map_.at(from) / mass_map_.at(to) * from_value;
  }
  double convert(const EnergyUnit& from, const EnergyUnit& to,
                 const double& from_value) const noexcept {
    return energy_map_.at(from) / energy_map_.at(to) * from_value;
  }
};
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_UNITCONVERTER_H
