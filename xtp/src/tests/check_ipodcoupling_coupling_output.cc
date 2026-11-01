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

// Checks the coupling output of ipodcoupling (ipodcoupling.jobs) against
// a reference.
//
// The H-saturation geometry depends on the OpenBabel version
// (RelaxNewAtoms; 3.1.1 vs. 3.2.0), and the DFT step amplifies these
// differences in the couplings. The run stage
// (xtp/src/tests/CMakeLists.txt) therefore starts from a stored orbFileAB
// with converged orbitals and runs only tasks="podcoupling", so the input
// is identical on all platforms. The geometry is checked separately in
// check_ipodcoupling_input_geometry.cc.
//
// Every coupling element must have a j attribute; j must agree with the
// reference up to an arbitrary sign (MO phase, cf. --signflip-attr j of
// votca_compare_xml) within a tolerance that only allows for rounding
// differences (compiler, BLAS, threading).

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE check_ipodcoupling_coupling_output_test
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>

// VOTCA includes
#include <votca/tools/property.h>

using namespace votca;

BOOST_AUTO_TEST_SUITE(check_ipodcoupling_coupling_output_test)

BOOST_AUTO_TEST_CASE(check_j_values_present_and_agree_with_reference) {

  // Paths are set at configure time (xtp/src/tests/CMakeLists.txt), since
  // command-line arguments would conflict with Boost.Test's argument
  // parsing (--run_test, --log_level, ...).
  std::string actual_file = IPODCOUPLING_RUN_JOBS_FILE;
  std::string reference_file = IPODCOUPLING_RUN_JOBS_REFERENCE_FILE;

  // Input is identical on all platforms, so only rounding differences
  // need to be tolerated.
  double tolerance_eV = 1e-4;

  tools::Property actual_xml;
  actual_xml.LoadFromXML(actual_file);
  std::vector<tools::Property *> actual_couplings =
      actual_xml.Select("jobs.job.output.ipodcoupling.coupling");

  tools::Property reference_xml;
  reference_xml.LoadFromXML(reference_file);
  std::vector<tools::Property *> reference_couplings =
      reference_xml.Select("jobs.job.output.ipodcoupling.coupling");

  // An empty reference would make all checks below meaningless.
  BOOST_REQUIRE_MESSAGE(!reference_couplings.empty(),
                        "Reference file '"
                            << reference_file
                            << "' contains no coupling elements.");

  // Keyed by (levelA, levelB), independent of the element order.
  std::map<std::pair<Index, Index>, double> reference_by_level;
  for (tools::Property *ref_coupling : reference_couplings) {
    Index levelA = ref_coupling->getAttribute<Index>("levelA");
    Index levelB = ref_coupling->getAttribute<Index>("levelB");
    double j_ref = ref_coupling->getAttribute<double>("j");
    reference_by_level[{levelA, levelB}] = j_ref;
  }

  BOOST_REQUIRE_MESSAGE(!actual_couplings.empty(),
                        "Output file '"
                            << actual_file
                            << "' contains no coupling elements; PODCoupling "
                               "produced no output.");

  for (tools::Property *actual_coupling : actual_couplings) {
    Index levelA = actual_coupling->getAttribute<Index>("levelA");
    Index levelB = actual_coupling->getAttribute<Index>("levelB");

    // A missing j means PODCoupling did not run or failed.
    BOOST_REQUIRE_MESSAGE(
        actual_coupling->hasAttribute("j"),
        "coupling element (levelA=" << levelA << ", levelB=" << levelB
                                    << ") has no j attribute.");

    double j_actual = actual_coupling->getAttribute<double>("j");

    auto it = reference_by_level.find({levelA, levelB});
    BOOST_REQUIRE_MESSAGE(
        it != reference_by_level.end(),
        "coupling element (levelA=" << levelA << ", levelB=" << levelB
                                    << ") has no matching reference entry.");

    double j_reference = it->second;
    // The sign of a coupling depends on the arbitrary MO phases; compare
    // with the smaller of the direct and the sign-flipped deviation.
    double deviation_eV = std::min(std::abs(j_actual - j_reference),
                                   std::abs(j_actual + j_reference));
    BOOST_CHECK_MESSAGE(
        deviation_eV <= tolerance_eV,
        "coupling (levelA=" << levelA << ", levelB=" << levelB
                            << ") deviates from the reference (allowing "
                               "for an arbitrary overall sign) by "
                            << deviation_eV << " eV (actual=" << j_actual
                            << ", reference=" << j_reference
                            << ", tolerance=" << tolerance_eV << " eV).");
  }
}

BOOST_AUTO_TEST_SUITE_END()
