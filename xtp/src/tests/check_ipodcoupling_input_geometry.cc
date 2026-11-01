/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Checks the H-saturated, OpenBabel-relaxed supermolecule geometry that
// the ipodcoupling "input" task writes to temp.inp. The test uses the
// ORCA package settings (ipodcoupling_orca_input.xml), since
// Orca::WriteInputFile writes a plain-text input file without running
// ORCA (XTPDFT keeps its input in memory). Together with
// check_ipodcoupling_coupling_output.cc (which starts from stored
// orbitals) this replaces a single, slow DFT+PODCoupling test.
//
// The geometry depends on the OpenBabel version (RelaxNewAtoms), but far
// less than couplings computed from it, which amplify geometry
// differences through the DFT calculation, so it can be checked with a
// small tolerance.

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE check_ipodcoupling_input_geometry_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Atom {
  std::string element;
  double x, y, z;
};

// Minimal parser for the Orca::WriteCoordinates format: a header line
// "* xyz charge spin", one line "element x y z" per atom (Angstrom), and
// a terminating "*" line. Whitespace and column widths are ignored.
std::vector<Atom> ParseOrcaInputGeometry(const std::string &filename) {
  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    throw std::runtime_error("Could not open '" + filename + "'");
  }

  std::vector<Atom> atoms;
  std::string line;
  bool in_geometry_block = false;
  while (std::getline(ifs, line)) {
    if (!in_geometry_block) {
      if (line.rfind("* xyz", 0) == 0) {
        in_geometry_block = true;
      }
      continue;
    }
    // End of the geometry block.
    std::string trimmed = line;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    if (!trimmed.empty() && trimmed[0] == '*') {
      break;
    }
    if (trimmed.empty()) {
      continue;
    }
    std::istringstream iss(line);
    Atom atom;
    if (iss >> atom.element >> atom.x >> atom.y >> atom.z) {
      atoms.push_back(atom);
    }
  }
  return atoms;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(check_ipodcoupling_input_geometry_test)

BOOST_AUTO_TEST_CASE(check_temp_inp_geometry_agrees_with_reference) {

  // Paths are set at configure time (xtp/src/tests/CMakeLists.txt), since
  // command-line arguments would conflict with Boost.Test's argument
  // parsing.
  std::string actual_file = IPODCOUPLING_INPUT_TEMP_INP_FILE;
  std::string reference_file = IPODCOUPLING_INPUT_TEMP_INP_REFERENCE_FILE;

  // Above the observed cross-platform displacement of up to about
  // 0.035 Angstrom (OpenBabel 3.2.0 on macOS vs. 3.1.1 on Ubuntu).
  double tolerance_angstrom = 0.05;

  std::vector<Atom> actual_atoms = ParseOrcaInputGeometry(actual_file);
  std::vector<Atom> reference_atoms = ParseOrcaInputGeometry(reference_file);

  BOOST_REQUIRE_MESSAGE(!reference_atoms.empty(),
                        "Reference file '"
                            << reference_file
                            << "' contains no parseable geometry.");
  BOOST_REQUIRE_MESSAGE(
      !actual_atoms.empty(),
      "File '" << actual_file
               << "' contains no parseable geometry; the ipodcoupling "
                  "\"input\" task did not write temp.inp.");

  BOOST_REQUIRE_MESSAGE(actual_atoms.size() == reference_atoms.size(),
                        "Atom count (" << actual_atoms.size()
                                       << ") does not match the reference atom "
                                          "count ("
                                       << reference_atoms.size() << ").");

  for (std::size_t i = 0; i < actual_atoms.size(); ++i) {
    const Atom &actual = actual_atoms[i];
    const Atom &reference = reference_atoms[i];

    BOOST_CHECK_MESSAGE(actual.element == reference.element,
                        "Atom " << i << ": element '" << actual.element
                                << "' does not match the reference element '"
                                << reference.element << "'.");

    double dx = actual.x - reference.x;
    double dy = actual.y - reference.y;
    double dz = actual.z - reference.z;
    double distance = std::sqrt(dx * dx + dy * dy + dz * dz);

    BOOST_CHECK_MESSAGE(
        distance <= tolerance_angstrom,
        "Atom " << i << " (" << actual.element
                << "): position deviates from the reference "
                   "by "
                << distance << " Angstrom (tolerance=" << tolerance_angstrom
                << " Angstrom) -- actual=(" << actual.x << ", " << actual.y
                << ", " << actual.z << "), reference=(" << reference.x << ", "
                << reference.y << ", " << reference.z << ").");
  }
}

BOOST_AUTO_TEST_SUITE_END()
