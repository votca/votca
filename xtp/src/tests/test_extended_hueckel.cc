/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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
#define BOOST_TEST_MODULE extended_hueckel_test

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "votca/xtp/extended_hueckel.h"

using namespace votca::xtp;

namespace {
constexpr double ev_to_ha = 0.03674932217565499;
}

BOOST_AUTO_TEST_SUITE(extended_hueckel_test)

BOOST_AUTO_TEST_CASE(has_and_has_element) {
  ExtendedHuckelParameters eht;

  BOOST_CHECK_EQUAL(eht.Has("H", 0), true);
  BOOST_CHECK_EQUAL(eht.Has("H", 1), false);

  BOOST_CHECK_EQUAL(eht.Has("C", 0), true);
  BOOST_CHECK_EQUAL(eht.Has("C", 1), true);
  BOOST_CHECK_EQUAL(eht.Has("C", 2), false);

  BOOST_CHECK_EQUAL(eht.Has("P", 0), true);
  BOOST_CHECK_EQUAL(eht.Has("P", 1), true);
  BOOST_CHECK_EQUAL(eht.Has("P", 2), true);

  BOOST_CHECK_EQUAL(eht.HasElement("H"), true);
  BOOST_CHECK_EQUAL(eht.HasElement("C"), true);
  BOOST_CHECK_EQUAL(eht.HasElement("Xe"), false);
}

BOOST_AUTO_TEST_CASE(get_exact_values) {
  ExtendedHuckelParameters eht;

  BOOST_CHECK_CLOSE(eht.Get("H", 0), -13.60 * ev_to_ha, 1e-10);
  BOOST_CHECK_CLOSE(eht.Get("C", 0), -21.40 * ev_to_ha, 1e-10);
  BOOST_CHECK_CLOSE(eht.Get("C", 1), -11.40 * ev_to_ha, 1e-10);
  BOOST_CHECK_CLOSE(eht.Get("P", 2), -14.00 * ev_to_ha, 1e-10);
  BOOST_CHECK_CLOSE(eht.Get("Se", 2), -14.80 * ev_to_ha, 1e-10);
}

BOOST_AUTO_TEST_CASE(get_throws_for_missing_exact_parameter) {
  ExtendedHuckelParameters eht;

  BOOST_REQUIRE_THROW(eht.Get("H", 1), std::runtime_error);
  BOOST_REQUIRE_THROW(eht.Get("C", 2), std::runtime_error);
  BOOST_REQUIRE_THROW(eht.Get("Xe", 0), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(get_with_fallback_exact_match) {
  ExtendedHuckelParameters eht;

  int used_l = -1;
  double eps = eht.GetWithFallback("C", 1, &used_l);

  BOOST_CHECK_CLOSE(eps, -11.40 * ev_to_ha, 1e-10);
  BOOST_CHECK_EQUAL(used_l, 1);
}

BOOST_AUTO_TEST_CASE(get_with_fallback_p_to_s) {
  ExtendedHuckelParameters eht;

  int used_l = -1;
  double eps = eht.GetWithFallback("H", 1, &used_l);

  BOOST_CHECK_CLOSE(eps, -13.60 * ev_to_ha, 1e-10);
  BOOST_CHECK_EQUAL(used_l, 0);
}

BOOST_AUTO_TEST_CASE(get_with_fallback_high_l_prefers_d_then_p_then_s) {
  ExtendedHuckelParameters eht;

  {
    int used_l = -1;
    double eps = eht.GetWithFallback("P", 5, &used_l);
    BOOST_CHECK_CLOSE(eps, -14.00 * ev_to_ha, 1e-10);
    BOOST_CHECK_EQUAL(used_l, 2);
  }

  {
    int used_l = -1;
    double eps = eht.GetWithFallback("C", 4, &used_l);
    BOOST_CHECK_CLOSE(eps, -11.40 * ev_to_ha, 1e-10);
    BOOST_CHECK_EQUAL(used_l, 1);
  }

  {
    int used_l = -1;
    double eps = eht.GetWithFallback("H", 3, &used_l);
    BOOST_CHECK_CLOSE(eps, -13.60 * ev_to_ha, 1e-10);
    BOOST_CHECK_EQUAL(used_l, 0);
  }
}

BOOST_AUTO_TEST_CASE(get_with_fallback_without_used_l_pointer) {
  ExtendedHuckelParameters eht;

  double eps = eht.GetWithFallback("Br", 4, nullptr);
  BOOST_CHECK_CLOSE(eps, -16.00 * ev_to_ha, 1e-10);
}

BOOST_AUTO_TEST_CASE(get_with_fallback_throws_for_missing_element) {
  ExtendedHuckelParameters eht;

  BOOST_REQUIRE_THROW(eht.GetWithFallback("Xe", 0), std::runtime_error);
  BOOST_REQUIRE_THROW(eht.GetWithFallback("U", 3), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()