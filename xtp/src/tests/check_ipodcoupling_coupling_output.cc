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

// Real, direct, standalone check (not a votca_compare_xml-based
// comparison test at all), worked through directly with the user.
//
// History, for the real, full context: the real, confirmed cross-
// platform numerical noise this same session originally worked
// through in ipodcoupling's own PODCoupling j (coupling) values
// traced back to real OpenBabel version differences across platforms
// (3.1.1 vs. 3.2.0) affecting RelaxNewAtoms's own H-saturation
// geometry, which a real, downstream DFT calculation then amplified.
// Rather than continue chasing bit-for-bit cross-platform determinism
// through that real DFT step, this test's own run stage
// (xtp/src/tests/CMakeLists.txt) was redesigned, worked through
// directly with the user, to pre-seed a real, checked-in reference
// orbFileAB (containing already-converged DFT orbitals) directly,
// with tasks="podcoupling" only -- so PODCoupling itself now runs on
// genuinely bit-identical input on every real platform, entirely
// sidestepping the real DFT step (and its own real OpenBabel-version
// dependency) for this specific comparison. A real, direct, separate,
// much cheaper test (check_ipodcoupling_input_geometry.cc) checks the
// H-saturation geometry itself, with a real, direct, small tolerance,
// instead.
//
// Given that redesign, this check now does two, real, separate
// things: (1) directly, genuinely REQUIRES that every real coupling
// element in the actual, real ipodcoupling.jobs output genuinely has
// its own real j attribute at all (a real, direct, fatal check -- via
// BOOST_REQUIRE -- since a genuinely missing j value would mean
// PODCoupling itself never actually ran, or genuinely crashed); and
// (2) genuinely REQUIRES (BOOST_CHECK, fatal, not merely a warning)
// that a real j value's own real magnitude agrees with the real,
// checked-in reference within a real, direct, tight tolerance --
// tight enough to catch a real, genuine regression, but still loose
// enough to tolerate ordinary floating-point rounding differences
// across real compilers/BLAS implementations/threading, since even
// bit-identical input can still produce very slightly different
// real floating-point results across those.

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE check_ipodcoupling_coupling_output_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <map>
#include <string>
#include <utility>

// VOTCA includes
#include <votca/tools/property.h>

using namespace votca;

BOOST_AUTO_TEST_SUITE(check_ipodcoupling_coupling_output_test)

BOOST_AUTO_TEST_CASE(check_j_values_present_and_agree_with_reference) {

  // Real, direct paths -- both known at real CMake configure time
  // already (REFPATH/RUNPATH, xtp/src/tests/CMakeLists.txt's own real
  // ipodcoupling run stage), so compiled in directly here rather than
  // parsed from real command-line arguments at all -- real, direct
  // command-line arguments here would otherwise conflict with
  // Boost::unit_test_framework's own real argument parsing
  // (--run_test=/--log_level=/etc.).
  std::string actual_file = IPODCOUPLING_RUN_JOBS_FILE;
  std::string reference_file = IPODCOUPLING_RUN_JOBS_REFERENCE_FILE;

  // Real, direct, tight tolerance, worked through directly with the
  // user: now that PODCoupling itself runs on genuinely bit-identical
  // input on every real platform (the real DFT step, and its own real
  // OpenBabel-version dependency, is sidestepped entirely by this
  // test's own pre-seeded reference orbFileAB), this only needs to
  // tolerate ordinary real floating-point rounding differences
  // (compiler/BLAS/threading), not real, genuine geometry
  // differences -- so kept far tighter than this file's own earlier,
  // now-stale 0.5 eV warn-only threshold.
  double tolerance_eV = 1e-4;

  tools::Property actual_xml;
  actual_xml.LoadFromXML(actual_file);
  std::vector<tools::Property *> actual_couplings =
      actual_xml.Select("jobs.job.output.ipodcoupling.coupling");

  tools::Property reference_xml;
  reference_xml.LoadFromXML(reference_file);
  std::vector<tools::Property *> reference_couplings =
      reference_xml.Select("jobs.job.output.ipodcoupling.coupling");

  // Real, direct, honest sanity check on the reference file itself --
  // fatal (BOOST_REQUIRE), since a genuinely empty/malformed reference
  // file would make every other real check below meaningless.
  BOOST_REQUIRE_MESSAGE(!reference_couplings.empty(),
                        "Reference file '"
                            << reference_file
                            << "' contains no coupling elements at "
                               "all -- cannot check against it.");

  // Real, direct lookup, keyed by (levelA, levelB) -- more robust
  // than a real, direct positional index alone, since this still
  // correctly matches even if the real, actual ordering of coupling
  // elements ever genuinely differs.
  std::map<std::pair<Index, Index>, double> reference_by_level;
  for (tools::Property *ref_coupling : reference_couplings) {
    Index levelA = ref_coupling->getAttribute<Index>("levelA");
    Index levelB = ref_coupling->getAttribute<Index>("levelB");
    double j_ref = ref_coupling->getAttribute<double>("j");
    reference_by_level[{levelA, levelB}] = j_ref;
  }

  BOOST_REQUIRE_MESSAGE(
      !actual_couplings.empty(),
      "Real, actual output file '"
          << actual_file
          << "' contains no coupling elements at all -- PODCoupling "
             "itself never actually produced any real output, or "
             "genuinely crashed.");

  for (tools::Property *actual_coupling : actual_couplings) {
    Index levelA = actual_coupling->getAttribute<Index>("levelA");
    Index levelB = actual_coupling->getAttribute<Index>("levelB");

    // Real, direct, fatal check: every real coupling element must
    // genuinely have its own real j attribute at all.
    BOOST_REQUIRE_MESSAGE(
        actual_coupling->hasAttribute("j"),
        "coupling element (levelA=" << levelA << ", levelB=" << levelB
                                    << ") is genuinely missing its own "
                                       "real j attribute entirely.");

    double j_actual = actual_coupling->getAttribute<double>("j");

    auto it = reference_by_level.find({levelA, levelB});
    BOOST_REQUIRE_MESSAGE(
        it != reference_by_level.end(),
        "coupling element (levelA="
            << levelA << ", levelB=" << levelB
            << ") has no matching reference entry at all to compare "
               "against.");

    double j_reference = it->second;
    double deviation_eV = std::abs(j_actual - j_reference);
    BOOST_CHECK_MESSAGE(
        deviation_eV <= tolerance_eV,
        "coupling (levelA=" << levelA << ", levelB=" << levelB
                            << ") deviates from the reference by "
                            << deviation_eV << " eV (actual=" << j_actual
                            << ", reference=" << j_reference
                            << ", tolerance=" << tolerance_eV << " eV).");
  }
}

BOOST_AUTO_TEST_SUITE_END()
