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
  grams_per_mole, // the same as atomic mass units
  kilograms,
  grams
};

enum TimeUnit { seconds, microseconds, nanoseconds, femtoseconds, picoseconds };

enum EnergyUnit {
  electron_volts,
  kilocalories,
  hartrees,
  joules,
  kilojoules,
  kilojoules_per_mole,
  joules_per_mole,
  kilocalories_per_mole
};

enum ChargeUnit { e, coulombs };

enum VelocityUnit {
  angstroms_per_femtosecond,
  angstroms_per_picosecond,
  nanometers_per_picosecond
};

enum ForceUnit {
  kilocalories_per_mole_angstrom,
  newtons,
  kilojoules_per_mole_nanometer,
  kilojoules_per_mole_angstrom,
  hatree_per_bohr
};
/**
 * @brief Class converts between different unit types
 */
class UnitConverter {
 private:
  /// All distances with respect to Ang
  constexpr double getDistanceValue_(const DistanceUnit& enum_type) const
      noexcept {
    switch (enum_type) {
      case DistanceUnit::meters:
        return 1E10;
      case DistanceUnit::centimeters:
        return 1E8;
      case DistanceUnit::nanometers:
        return 10.0;
      case DistanceUnit::angstroms:
        return 1.0;
      case DistanceUnit::bohr:
        return 0.52917721092;
    }
    return 0.0;
  }

  /// All times with respect to pico seconds
  constexpr double getTimeValue_(const TimeUnit& enum_type) const noexcept {
    switch (enum_type) {
      case TimeUnit::seconds:
        return 1E12;
      case TimeUnit::microseconds:
        return 1E6;
      case TimeUnit::nanoseconds:
        return 1000;
      case TimeUnit::picoseconds:
        return 1.0;
      case TimeUnit::femtoseconds:
        return 0.001;
    }
    return 0.0;
  }

  /// All masses with respect to atomic mass units
  constexpr double getMassValue_(const MassUnit& enum_type) const noexcept {
    switch (enum_type) {
      case MassUnit::kilograms:
        return 6.0221366517E26;
      case MassUnit::grams:
        return 6.0221366517E23;
      case MassUnit::picograms:
        return 6.0221366517E14;
      case MassUnit::femtograms:
        return 6.0221366517E11;
      case MassUnit::attograms:
        return 60221366517;
      case MassUnit::atomic_mass_units:
        return 1.0;
      case MassUnit::grams_per_mole:
        return 1.0;
    }
    return 0.0;
  }

  /// All energies in terms of electron volts
  constexpr double getEnergyValue_(const EnergyUnit& enum_type) const noexcept {
    switch (enum_type) {
      case EnergyUnit::kilojoules_per_mole:
        return 96.4853074993;
      case EnergyUnit::joules_per_mole:
        return 96.4853074993E3;
      case EnergyUnit::kilocalories_per_mole:
        return 23.061;
      case EnergyUnit::kilocalories:
        return 2.613195131836172E22;
      case EnergyUnit::joules:
        return 6.2415097E18;
      case EnergyUnit::hartrees:
        return 27.211368602;
      case EnergyUnit::electron_volts:
        return 1.0;
    }
    return 0.0;
  }

  /// All charge in terms of elementary charge e
  constexpr double getChargeValue_(const ChargeUnit& enum_type) const noexcept {
    switch (enum_type) {
      case ChargeUnit::e:
        return 1.0;
      case ChargeUnit::coulombs:
        return 1.602176565E-19;
    }
    return 0.0;
  }

 public:
  constexpr double convert(const DistanceUnit& from,
                           const DistanceUnit& to) const noexcept {
    return getDistanceValue_(from) / getDistanceValue_(to);
  }
  constexpr double convert(const TimeUnit& from, const TimeUnit& to) const
      noexcept {
    return getTimeValue_(from) / getTimeValue_(to);
  }
  constexpr double convert(const MassUnit& from, const MassUnit& to) const
      noexcept {
    return getMassValue_(from) / getMassValue_(to);
  }
  constexpr double convert(const EnergyUnit& from, const EnergyUnit& to) const
      noexcept {
    return getEnergyValue_(from) / getEnergyValue_(to);
  }
};
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_UNITCONVERTER_H
