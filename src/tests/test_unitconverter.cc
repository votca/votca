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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE unitconverter_test
#include <iostream>
#include "../../include/votca/tools/unitconverter.h"
#include <boost/test/unit_test.hpp>
using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(unitconverter_test)

BOOST_AUTO_TEST_CASE(unitconverter_test_distance) {
  UnitConverter converter;

  double distance = 1;  // Assuming in meters
  double distance_new =
      converter.convert(DistanceUnit::meters, DistanceUnit::nanometers) *
      distance;
  BOOST_CHECK_EQUAL(1E9, distance_new);

  distance_new =
      converter.convert(DistanceUnit::meters, DistanceUnit::centimeters) *
      distance;
  BOOST_CHECK_EQUAL(100, distance_new);

  distance_new =
      converter.convert(DistanceUnit::meters, DistanceUnit::angstroms) *
      distance;
  BOOST_CHECK_EQUAL(1E10, distance_new);

  distance_new =
      converter.convert(DistanceUnit::meters, DistanceUnit::bohr) * distance;
  BOOST_CHECK_CLOSE(1.8897e+10, distance_new, 0.01);

  // Now assuming distance in angstroms
  distance_new =
      converter.convert(DistanceUnit::angstroms, DistanceUnit::nanometers) * distance;
  BOOST_CHECK_CLOSE(0.1, distance_new, 0.01);

  distance_new =
      converter.convert(DistanceUnit::angstroms, DistanceUnit::meters) * distance;
  BOOST_CHECK_CLOSE(1.E-10, distance_new, 0.01);

  // Now assuming distance in Bohr
  distance_new =
      converter.convert(DistanceUnit::nanometers, DistanceUnit::bohr) * distance;
  BOOST_CHECK_CLOSE(18.8973, distance_new, 0.01);

}

BOOST_AUTO_TEST_CASE(unitconverter_test_time) {
  UnitConverter converter;

  double time = 1;  // Assuming in seconds
  double time_new =
      converter.convert(TimeUnit::seconds, TimeUnit::microseconds) * time;
  BOOST_CHECK_CLOSE(1E6, time_new, 0.01);

  time_new = converter.convert(TimeUnit::seconds, TimeUnit::nanoseconds) * time;
  BOOST_CHECK_CLOSE(1E9, time_new, 0.01);

  time_new = converter.convert(TimeUnit::seconds, TimeUnit::picoseconds) * time;
  BOOST_CHECK_CLOSE(1E12, time_new, 0.01);

  time_new =
      converter.convert(TimeUnit::seconds, TimeUnit::femtoseconds) * time;
  BOOST_CHECK_CLOSE(1E15, time_new, 0.01);
}

BOOST_AUTO_TEST_CASE(unitconverter_test_mass) {
  UnitConverter converter;

  double mass = 1;  // Assuming in kilograms
  double mass_new =
      converter.convert(MassUnit::kilograms, MassUnit::grams) * mass;
  BOOST_CHECK_CLOSE(1E3, mass_new, 0.01);

  mass_new = converter.convert(MassUnit::kilograms, MassUnit::picograms) * mass;
  BOOST_CHECK_CLOSE(1E15, mass_new, 0.01);

  mass_new =
      converter.convert(MassUnit::kilograms, MassUnit::femtograms) * mass;
  BOOST_CHECK_CLOSE(1E18, mass_new, 0.01);

  mass_new = converter.convert(MassUnit::kilograms, MassUnit::attograms) * mass;
  BOOST_CHECK_CLOSE(1E21, mass_new, 0.01);

  mass_new =
      converter.convert(MassUnit::kilograms, MassUnit::atomic_mass_units) *
      mass;
  BOOST_CHECK_CLOSE(6.02213665E26, mass_new, 0.01);

  mass_new =
      converter.convert(MassUnit::kilograms, MassUnit::grams_per_mole) *
      mass;
  BOOST_CHECK_CLOSE(6.02213665E26, mass_new, 0.01);
}

BOOST_AUTO_TEST_CASE(unitconverter_test_energy) {
  UnitConverter converter;

  double energy = 1;  // Assuming in electron volts
  double energy_new = converter.convert(EnergyUnit::electron_volts,
                                        EnergyUnit::electron_volts) *
                      energy;
  BOOST_CHECK_CLOSE(1, energy_new, 0.01);

  energy_new =
      converter.convert(EnergyUnit::electron_volts, EnergyUnit::kilocalories) *
      energy;
  BOOST_CHECK_CLOSE(3.82929E-23, energy_new, 0.01);

  energy_new =
      converter.convert(EnergyUnit::electron_volts, EnergyUnit::joules) *
      energy;
  BOOST_CHECK_CLOSE(1.60218E-19, energy_new, 0.01);

  energy_new =
      converter.convert(EnergyUnit::electron_volts, EnergyUnit::hartrees) *
      energy;
  BOOST_CHECK_CLOSE(0.0367493, energy_new, 0.01);

  energy_new =
      converter.convert(EnergyUnit::electron_volts, EnergyUnit::kilojoules) *
      energy;
  BOOST_CHECK_CLOSE(1.602176E-22, energy_new, 0.01);
}

BOOST_AUTO_TEST_CASE(unitconverter_test_molar_energy) {
  UnitConverter converter;

  double energy = 1;  // Assuming in electron volts
  double energy_new = converter.convert(MolarEnergyUnit::electron_volts_per_mole,
                                        MolarEnergyUnit::electron_volts_per_mole) *
                      energy;
  BOOST_CHECK_CLOSE(1, energy_new, 0.01);

  energy_new =
      converter.convert(MolarEnergyUnit::electron_volts_per_mole, MolarEnergyUnit::kilocalories_per_mole) *
      energy;
  BOOST_CHECK_CLOSE(3.82929E-23, energy_new, 0.01);

  energy_new =
      converter.convert(MolarEnergyUnit::electron_volts_per_mole, MolarEnergyUnit::joules_per_mole) *
      energy;
  BOOST_CHECK_CLOSE(1.60218E-19, energy_new, 0.01);

  energy_new =
      converter.convert(MolarEnergyUnit::electron_volts_per_mole, MolarEnergyUnit::hartrees_per_mole) *
      energy;
  BOOST_CHECK_CLOSE(0.0367493, energy_new, 0.01);

  energy_new =
      converter.convert(MolarEnergyUnit::electron_volts_per_mole, MolarEnergyUnit::kilojoules_per_mole) *
      energy;
  BOOST_CHECK_CLOSE(1.602176E-22, energy_new, 0.01);
}

BOOST_AUTO_TEST_CASE(unitconverter_test_charge) {
  UnitConverter converter;

  double charge = 1.602176565E-19;  // Assuming in units of coulombs
  double charge_new = converter.convert(ChargeUnit::coulombs,
                                        ChargeUnit::e) * charge;
  BOOST_CHECK_CLOSE(1.000, charge_new, 0.01 );

  charge = 1; // Assuming in eV
  charge_new = converter.convert(ChargeUnit::e,ChargeUnit::coulombs) * charge;
  BOOST_CHECK_CLOSE(1.60217656E-19, charge_new, 0.01);
}

BOOST_AUTO_TEST_CASE(unitconverter_test_velocity) {
  UnitConverter converter;

  double velocity = 1;  // Assuming in units of angstroms_per_femtosecond
  double velocity_new = converter.convert(VelocityUnit::angstroms_per_femtosecond,
                                        VelocityUnit::angstroms_per_picosecond) * velocity;
  BOOST_CHECK_CLOSE(1000, velocity_new, 0.01 );

  velocity_new = converter.convert(VelocityUnit::angstroms_per_femtosecond,
      VelocityUnit::nanometers_per_picosecond) * velocity;
  BOOST_CHECK_CLOSE(100, velocity_new, 0.01 );

  // Assuming in units of nanometers_per_picosecond
  velocity_new = converter.convert(VelocityUnit::nanometers_per_picosecond,
      VelocityUnit::angstroms_per_femtosecond) * velocity;
  BOOST_CHECK_CLOSE(0.01, velocity_new, 0.01 );
  
}

BOOST_AUTO_TEST_CASE(unitconverter_test_force) {
  UnitConverter converter;

  double force = 1;  // Assuming kilojoules_per_nanometer
  double force_new = converter.convert(ForceUnit::kilojoules_per_nanometer,
                                        ForceUnit::kilocalories_per_angstrom) * force;
  BOOST_CHECK_CLOSE(0.0239006, force_new, 0.01 );

  force_new = converter.convert(ForceUnit::kilojoules_per_nanometer,
      ForceUnit::kilojoules_per_angstrom) * force;
  BOOST_CHECK_CLOSE(0.1, force_new, 0.01 );

  force_new = converter.convert(ForceUnit::kilojoules_per_nanometer, ForceUnit::hatree_per_bohr) * force;
  BOOST_CHECK_CLOSE(1.2138E19, force_new, 0.01 );

  force_new = converter.convert(ForceUnit::kilojoules_per_nanometer,
      ForceUnit::newtons) * force;
  BOOST_CHECK_CLOSE(1E12, force_new, 0.01 );
}

BOOST_AUTO_TEST_CASE(unitconverter_test_molar_force) {
  UnitConverter converter;

  double force = 1;  // Assuming kilojoules_per_mole_nanometer
  double force_new = converter.convert(MolarForceUnit::kilojoules_per_mole_nanometer,
                                        MolarForceUnit::kilocalories_per_mole_angstrom) * force;
  BOOST_CHECK_CLOSE(0.0239006, force_new, 0.01 );

  force_new = converter.convert(MolarForceUnit::kilojoules_per_mole_nanometer,
      MolarForceUnit::kilojoules_per_mole_angstrom) * force;
  BOOST_CHECK_CLOSE(0.1, force_new, 0.01 );

  force_new = converter.convert(MolarForceUnit::kilojoules_per_mole_nanometer, MolarForceUnit::hatree_per_mole_bohr) * force;
  BOOST_CHECK_CLOSE(1.2138E19, force_new, 0.01 );

  force_new = converter.convert(MolarForceUnit::kilojoules_per_mole_nanometer,
      MolarForceUnit::newtons_per_mole) * force;
  BOOST_CHECK_CLOSE(1E12, force_new, 0.01 );
}
BOOST_AUTO_TEST_SUITE_END()
