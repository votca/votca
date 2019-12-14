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

// Extrinsic Energy Units
enum EnergyUnit {
  electron_volts,
  kilocalories,
  hartrees,
  joules,
  kilojoules,
};

// Intrinsic Energy Units
enum MolarEnergyUnit {
  kilojoules_per_mole,
  joules_per_mole,
  kilocalories_per_mole,
  electron_volts_per_mole,
  hartrees_per_mole
};

enum ChargeUnit { e, coulombs };

enum VelocityUnit {
  angstroms_per_femtosecond,
  angstroms_per_picosecond,
  nanometers_per_picosecond
};

// Extrinsic Force Units
enum ForceUnit {
  kilocalories_per_angstrom,
  newtons,
  kilojoules_per_nanometer,
  kilojoules_per_angstrom,
  hatree_per_bohr
};

// Intrinsic Force Units
enum MolarForceUnit {
  kilocalories_per_mole_angstrom,
  newtons_per_mole,
  kilojoules_per_mole_nanometer,
  kilojoules_per_mole_angstrom,
  hatree_per_mole_bohr
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
        return 1E-10;
      case DistanceUnit::centimeters:
        return 1E-8;
      case DistanceUnit::nanometers:
        return 0.1;
      case DistanceUnit::angstroms:
        return 1.0;
      case DistanceUnit::bohr:
        return 1.8897161646321;
    }
    return 0.0;
  }

  /// All times with respect to pico seconds
  constexpr double getTimeValue_(const TimeUnit& enum_type) const noexcept {
    switch (enum_type) {
      case TimeUnit::seconds:
        return 1E-12;
      case TimeUnit::microseconds:
        return 1E-6;
      case TimeUnit::nanoseconds:
        return 1E-3;
      case TimeUnit::picoseconds:
        return 1.0;
      case TimeUnit::femtoseconds:
        return 1000;
    }
    return 0.0;
  }

  /// All masses with respect to atomic mass units
  constexpr double getMassValue_(const MassUnit& enum_type) const noexcept {
    switch (enum_type) {
      case MassUnit::kilograms:
        return 1.6605402E-27;
      case MassUnit::grams:
        return 1.6605402E-24;
      case MassUnit::picograms:
        return 1.6605402E-12;
      case MassUnit::femtograms:
        return 1.6605402E-9;
      case MassUnit::attograms:
        return 1.66054019E-6;
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
      case EnergyUnit::kilocalories:
        return 3.82929389E-23;
      case EnergyUnit::kilojoules:
        return 1.602176634E-22;
      case EnergyUnit::joules:
        return 1.602176634E-19;
      case EnergyUnit::hartrees:
        return 0.0367493;
      case EnergyUnit::electron_volts:
        return 1.0;
    }
    return 0.0;
  }

  // All molar energies in terms of eV_per_mole
  constexpr double getMolarEnergyValue_(const MolarEnergyUnit& enum_type) const noexcept {
    switch (enum_type) {
      case MolarEnergyUnit::kilojoules_per_mole:
        return 1.602176634E-22;
      case MolarEnergyUnit::joules_per_mole:
        return 1.602176634E-19;
      case MolarEnergyUnit::kilocalories_per_mole:
        return 3.82929389E-23;
      case MolarEnergyUnit::hartrees_per_mole:
        return 0.0367493;
      case MolarEnergyUnit::electron_volts_per_mole:
        return 1.0;
    }
    return 0.0;
  }
  /// All charge in terms of elementary charge e
  constexpr double getChargeValue_(const ChargeUnit& enum_type) const noexcept {
    switch (enum_type) {
      case ChargeUnit::e:
        return 1;
      case ChargeUnit::coulombs:
        return 1.602176565E-19;
    }
    return 0.0;
  }

  /// All velocities in terms of nanometers per picosecond
  constexpr double getVelocityValue_(const VelocityUnit& enum_type) const
    noexcept {
      switch (enum_type) {
        case VelocityUnit::nanometers_per_picosecond:
          return 1.0;
        case VelocityUnit::angstroms_per_picosecond:
          return convert(DistanceUnit::nanometers, DistanceUnit::angstroms);
        case VelocityUnit::angstroms_per_femtosecond:
          return convert(DistanceUnit::nanometers, DistanceUnit::angstroms) / convert(TimeUnit::picoseconds, TimeUnit::femtoseconds);
      }
      return 0.0;
    }

   /// Default force unit is the kilojoules per nanometer
  constexpr double getForceValue_(const ForceUnit& enum_type) const noexcept {
    switch (enum_type) {
      case ForceUnit::kilocalories_per_angstrom:
        return convert(EnergyUnit::kilojoules,
            EnergyUnit::kilocalories) /
          convert(DistanceUnit::nanometers, DistanceUnit::angstroms);
      case ForceUnit::newtons:
        return convert(EnergyUnit::kilojoules, EnergyUnit::joules) /
          convert(DistanceUnit::nanometers, DistanceUnit::meters);
      case ForceUnit::kilojoules_per_nanometer:
        return 1.0;
      case ForceUnit::kilojoules_per_angstrom:
        return (1.0) / convert(DistanceUnit::nanometers, DistanceUnit::angstroms);
      case ForceUnit::hatree_per_bohr:
        return convert(EnergyUnit::kilojoules, EnergyUnit::hartrees) /
        (convert(DistanceUnit::nanometers, DistanceUnit::bohr));
    }
    return 0.0;
  }

   /// Default force unit is the kilojoules per mole nanometer
  constexpr double getMolarForceValue_(const MolarForceUnit& enum_type) const noexcept {
    switch (enum_type) {
      case MolarForceUnit::kilocalories_per_mole_angstrom:
        return convert(MolarEnergyUnit::kilojoules_per_mole,
            MolarEnergyUnit::kilocalories_per_mole) /
          convert(DistanceUnit::nanometers, DistanceUnit::angstroms);
      case MolarForceUnit::newtons_per_mole:
        return convert(MolarEnergyUnit::kilojoules_per_mole, MolarEnergyUnit::joules_per_mole) /
          convert(DistanceUnit::nanometers, DistanceUnit::meters);
      case MolarForceUnit::kilojoules_per_mole_nanometer:
        return 1.0;
      case MolarForceUnit::kilojoules_per_mole_angstrom:
        return (1.0) / convert(DistanceUnit::nanometers, DistanceUnit::angstroms);
      case MolarForceUnit::hatree_per_mole_bohr:
        return convert(MolarEnergyUnit::kilojoules_per_mole, MolarEnergyUnit::hartrees_per_mole) /
        (convert(DistanceUnit::nanometers, DistanceUnit::bohr));
    }
    return 0.0;
  }
 public:
  constexpr double convert(const DistanceUnit& from,
                           const DistanceUnit& to) const noexcept {
    return getDistanceValue_(to) / getDistanceValue_(from);
  }
  constexpr double convert(const TimeUnit& from, const TimeUnit& to) const
      noexcept {
    return getTimeValue_(to) / getTimeValue_(from);
  }
  constexpr double convert(const MassUnit& from, const MassUnit& to) const
      noexcept {
    return getMassValue_(to) / getMassValue_(from);
  }
  constexpr double convert(const EnergyUnit& from, const EnergyUnit& to) const noexcept {
    return getEnergyValue_(to) / getEnergyValue_(from);
  }
  constexpr double convert(const MolarEnergyUnit& from, const MolarEnergyUnit& to) const noexcept {
    return getMolarEnergyValue_(to) / getMolarEnergyValue_(from);
  }
  constexpr double convert(const ChargeUnit& from, const ChargeUnit& to) const
    noexcept {
      return getChargeValue_(to) / getChargeValue_(from);
    }
  constexpr double convert(const VelocityUnit& from,
      const VelocityUnit& to) const noexcept {
    return getVelocityValue_(to) / getVelocityValue_(from);
  }
  constexpr double convert(const ForceUnit& from, const ForceUnit& to) const
    noexcept {
      return (getForceValue_(to)) / (getForceValue_(from));
    }
  constexpr double convert(const MolarForceUnit& from, const MolarForceUnit& to) const
    noexcept {
      return (getMolarForceValue_(to)) / (getMolarForceValue_(from));
    }

};
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_UNITCONVERTER_H
