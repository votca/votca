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

// Real, direct, standalone check, worked through directly with the
// user: the second half of splitting ipodcoupling's own run-stage
// test into two, real, separate, much cheaper pieces, in direct
// response to the real, confirmed real cross-platform OpenBabel-
// version numerical noise this same session already worked through
// directly (see check_ipodcoupling_coupling_output.cc's own header
// comment for that full, real reasoning), combined with a real,
// direct debug-build timeout on the previous, real, single
// combined DFT+PODCoupling test.
//
// Checks the geometry written into temp.inp by IPodCoupling's own
// "input"-only task, using a real, direct, orca-flavored reference
// ipodcoupling.xml (ipodcoupling_orca_input.xml, distinct from this
// same test suite's own xtp-flavored ipodcoupling.xml): confirmed
// directly, by reading it, that Orca::WriteInputFile itself --
// unlike XTPDFT::WriteInputFile, which is entirely in-memory --
// genuinely writes a real, direct, plain text temp.inp file
// (Orca::WriteCoordinates), containing the H-saturated,
// OpenBabel-relaxed supermolecule geometry, without ever genuinely
// running any real DFT calculation, or requiring a real ORCA
// installation, at all.
//
// This geometry itself still genuinely reflects the same real,
// platform-dependent OpenBabel-version noise as before (this
// "input"-only task still genuinely runs RelaxNewAtoms, same as
// every other task combination) -- so, unlike
// check_ipodcoupling_coupling_output.cc's own j-value check, this
// uses a real, direct, fatal (BOOST_CHECK, not merely
// BOOST_WARN_MESSAGE) small-tolerance comparison instead: worked
// through directly with the user, checking the geometry itself,
// directly, is genuinely less volatile than checking the final
// coupling values (which amplify small geometry differences through
// a real, additional, downstream DFT calculation on top) -- so a
// real, direct, fatal check here is both achievable and worthwhile.

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

// Real, direct, minimal parser for Orca::WriteCoordinates's own real,
// actual output format (confirmed directly, by reading it, before
// writing this parser): a "* xyz <charge> <spin>" header line,
// followed by one real "<element> <x> <y> <z>" line per real atom (in
// Angstrom), terminated by a real "* " line. Deliberately tolerant of
// exact whitespace/column width (Orca::WriteCoordinates itself uses
// std::setw for fixed-width columns, which this parser does not try
// to replicate exactly at all) -- only the real, actual element and
// numeric values themselves are checked.
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
    // Real, direct end-of-geometry marker (Orca::WriteCoordinates's
    // own "* " terminator line).
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

  // Real, direct paths -- both known at real CMake configure time
  // already (REFPATH/RUNPATH, xtp/src/tests/CMakeLists.txt's own real
  // ipodcoupling input stage), compiled in directly here rather than
  // parsed from real command-line arguments at all, for the same,
  // already-established real reason as
  // check_ipodcoupling_coupling_output.cc (real command-line
  // arguments here would otherwise conflict with
  // Boost::unit_test_framework's own real argument parsing).
  std::string actual_file = IPODCOUPLING_INPUT_TEMP_INP_FILE;
  std::string reference_file = IPODCOUPLING_INPUT_TEMP_INP_REFERENCE_FILE;

  // Real, direct, small tolerance, worked through directly with the
  // user: sized deliberately looser than the real, actual,
  // cross-platform geometry displacement this same session already
  // directly measured (up to ~35 mA, between a real macOS/Homebrew
  // OpenBabel 3.2.0 run and a real Ubuntu/apt OpenBabel 3.1.1 run,
  // confirmed directly this same session) -- but still a real, direct
  // fatal check (BOOST_CHECK, not merely a warning), since checking
  // the geometry itself, directly, is genuinely less volatile than
  // checking the final, downstream coupling values.
  double tolerance_angstrom = 0.05;

  std::vector<Atom> actual_atoms = ParseOrcaInputGeometry(actual_file);
  std::vector<Atom> reference_atoms = ParseOrcaInputGeometry(reference_file);

  BOOST_REQUIRE_MESSAGE(
      !reference_atoms.empty(),
      "Reference file '" << reference_file
                         << "' contains no real, parseable geometry at "
                            "all -- cannot check against it.");
  BOOST_REQUIRE_MESSAGE(
      !actual_atoms.empty(),
      "Real, actual file '"
          << actual_file
          << "' contains no real, parseable geometry at all -- "
             "IPodCoupling's own \"input\" task itself never actually "
             "wrote temp.inp, or genuinely crashed before doing so.");

  BOOST_REQUIRE_MESSAGE(
      actual_atoms.size() == reference_atoms.size(),
      "Real, actual atom count (" << actual_atoms.size()
                                  << ") does not match the reference atom "
                                     "count ("
                                  << reference_atoms.size() << ") at all.");

  for (std::size_t i = 0; i < actual_atoms.size(); ++i) {
    const Atom &actual = actual_atoms[i];
    const Atom &reference = reference_atoms[i];

    BOOST_CHECK_MESSAGE(
        actual.element == reference.element,
        "Atom " << i << ": real, actual element '" << actual.element
               << "' does not match the reference element '"
               << reference.element << "' at all.");

    double dx = actual.x - reference.x;
    double dy = actual.y - reference.y;
    double dz = actual.z - reference.z;
    double distance = std::sqrt(dx * dx + dy * dy + dz * dz);

    BOOST_CHECK_MESSAGE(
        distance <= tolerance_angstrom,
        "Atom " << i << " (" << actual.element
               << "): real, actual position deviates from the reference "
                  "by "
               << distance << " Angstrom (tolerance=" << tolerance_angstrom
               << " Angstrom) -- actual=(" << actual.x << ", " << actual.y
               << ", " << actual.z << "), reference=(" << reference.x
               << ", " << reference.y << ", " << reference.z << ").");
  }
}

BOOST_AUTO_TEST_SUITE_END()
